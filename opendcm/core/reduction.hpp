/*
    openDCM, dimensional constraint manager
    Copyright (C) 2014  Stefan Troeger <stefantroeger@gmx.net>

    This library is free software; you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License along
    with this library; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/

#ifndef DCM_ANALYSE_H
#define DCM_ANALYSE_H

#include <Eigen/Core>

#include "geometry.hpp"
#include "clustergraph.hpp"
#include "constraint.hpp"
#include "scheduler.hpp"
#include "utilities.hpp"

#include <boost/multi_array.hpp>
#include <unordered_map>
#include <typeindex>
#include <vector>
                

namespace dcm {
namespace numeric {   

/**
 * @brief Structure to setup the numeric equations for a solver
 *
 * This struct is the entry point for the numeric solver setup. It has all the functions needed to get the
 * numeric geometries and constraints.
 */
template<typename Kernel>
struct EquationBuilder {
    
    virtual ~EquationBuilder(){};
    
    typedef shedule::FlowGraph::Node  Node;
    typedef std::shared_ptr<numeric::Calculatable<Kernel>> CalcPtr;
    
    //this function is used to create a default numeric geometry for the given symbolic one.
    virtual CalcPtr createGeometry(const graph::LocalVertex& vertex) = 0;
     
    //this function creates a default numeric geometry equation and a node for it in the 
    //corresponding flow graph
    virtual std::pair< CalcPtr, Node>
    createGeometryNode(const graph::LocalVertex& vertex, shedule::FlowGraph& flowgraph) = 0;
    
    //this function creates a reduced numeric geometry equation. For this the equation it depends on 
    //is required
    virtual CalcPtr createReducedGeometry(const graph::LocalVertex& vertex, CalcPtr basegeometry) = 0;
     
    //this function is used to create a reduced numeric geometry equation and a node for it in the 
    //corresponding flow graph. For this the equation it depends on is required
    virtual std::pair< CalcPtr, Node>
    createReducedGeometryNode(const graph::LocalVertex& vertex, CalcPtr basegeometry,
                              shedule::FlowGraph& flowgraph) = 0;
    
    //this function is used to create all default constraint equations
    virtual std::vector<CalcPtr>
    createEquations(CalcPtr g1, CalcPtr g2) = 0;
                       
    //this function is used to create a node in the calculation flow graph with all remaining
    //constraint Equations
    virtual std::pair<std::vector<CalcPtr>, Node>
    createEquationsNode(CalcPtr g1, CalcPtr g2, shedule::FlowGraph& flow) = 0;

    //this function is used to create all default constraint equations
    virtual std::vector<CalcPtr>
    createReducedEquations(CalcPtr g1, CalcPtr g2) = 0;
                       
    //this function is used to create a node in the calculation flow graph with all remaining
    //constraint Equations
    virtual std::pair<std::vector<CalcPtr>, Node>
    createReducedEquationsNode(CalcPtr g1, CalcPtr g2, shedule::FlowGraph& flow) = 0;
};

template<typename Kernel>
struct EquationBuilderProperty {
    typedef EquationBuilder<Kernel>* type;
    struct default_value {
        EquationBuilder<Kernel>* operator()() {
            return nullptr;
        };
    };
};

}//numeric

/**
 * @brief The tree structure used for Constrain reduction
 * 
 * 
 */
namespace reduction {
   
struct Edge;
struct Node;

/**
 * @brief Data structure for tree traversal
 * 
 * As the reduction tree should be traversed in parallel the nodes and edges must be reentrant, hence
 * they do not hold any data and only provide algorithms. The data for the traversal is stored in the
 * treewalker.
 */
struct TreeWalker {
    
    Node* getInitialNode()      {return m_initialNode;};
    void  setInitialNode(Node* n) {m_initialNode = n;};

    Node* getFinalNode()      {return m_finalNode;};
    void  setFinalNode(Node* n) {m_finalNode = n;};
    
private:
    Node* m_initialNode;
    Node* m_finalNode;
};

template<typename Kernel>
struct ConstraintWalker : public TreeWalker {
    
    typedef std::shared_ptr<numeric::Calculatable<Kernel>>                Equation;
    typedef typename std::vector<symbolic::Constraint*>::const_iterator   iterator;  
    
    void setConstraintPool(std::vector<symbolic::Constraint*> c) {m_constraintPool = c;};
    
    //interface for accessing constraints
    bool     empty() const {return m_constraintPool.empty();}
    int      size()  const {return m_constraintPool.size();}
    iterator begin() const {return m_constraintPool.begin();};
    iterator end()   const {return m_constraintPool.end();};
    
    //multiple constraints of the same type are always reduntant or conflicting, they do not occure here
    template<typename T>
    symbolic::TypeConstraint<T>* getConstraint(int type) {
        for(auto c : m_constraintPool) {
            if(c->type == type)
                return static_cast<symbolic::TypeConstraint<T>*>(c);
        }
        return nullptr;
    };
    void acceptConstraint(symbolic::Constraint* c) {
        m_constraintPool.erase(std::remove(m_constraintPool.begin(), m_constraintPool.end(), c), m_constraintPool.end());
    };
    
private:
    std::vector<symbolic::Constraint*>  m_constraintPool;  //all the constraints we want to reduce
};

template<typename Kernel, template<class> class Primitive>
struct TargetWalker : public ConstraintWalker<Kernel> {
    
    typedef ConstraintWalker<Kernel>                                      Inherited;
    typedef std::shared_ptr<numeric::Equation<Kernel, Primitive<Kernel>>> Geometry;
    typedef std::shared_ptr<numeric::Calculatable<Kernel>>                Equation;
    
    TargetWalker(const Primitive<Kernel>& prim) : m_targetPrimitive(prim) {};
 
    const Primitive<Kernel>& getTargetPrimitive() {return m_targetPrimitive;};
        
    void        setInputEquation(Equation e) {m_inputEqn = e;};
    Equation    getInputEquation() {return m_inputEqn;};
    
private:
    const Primitive<Kernel>&    m_targetPrimitive; //primitive holding the value
    Equation                    m_inputEqn;  //cummulative input equation
};

template<typename Kernel, template<class> class SourcePrimitive, template<class> class TargetPrimitive>
struct SourceTargetWalker : TargetWalker<Kernel, TargetPrimitive>{
  
    SourceTargetWalker(const SourcePrimitive<Kernel>& sprim, const TargetPrimitive<Kernel>& tprim) 
                            : TargetWalker<Kernel, TargetPrimitive>(tprim), m_sourcePrimitive(sprim) {};
    
    const SourcePrimitive<Kernel>& getSourcePrimitive() {return m_sourcePrimitive;};

private:
    const SourcePrimitive<Kernel>&    m_sourcePrimitive; //primitive holding the value

};

/**
 * @brief Connects two nodes in a conditional manner
 *
 * This class describes the connection between two Nodes in the reduction tree. Edges are conditional,
 * they are only used as connections if some criteria is meet. To check this criteria the method 
 * \ref apply is used, it returns the validity of the connection. Furthermore an edge can have a own
 * behavior, an action can be executed. For this is custom class has to be derived which overrides 
 * the provided apply function.
 */
struct Edge {

    Node* source;
    Node* target;

    virtual ~Edge() {};
    
    /**
     * @brief Checks edge condition and traverses the tree
     * 
     * This class handles all 3 tasks of the edge: It checks for validity and returns true/false 
     * dependent on the criteria, it executes any user defined action when overridden in a derived
     * class and it further executes the tree traversal. Note that this base version has no criteria
     * assigned, it always returns true. It executes no action and only further traverses the tree.
     * \Note When overriding this function it must be assured that the base version is called to ensure
     * correct tree traversing.
     * 
     * @param walker The TreeWalker holding all data
     * @return bool  True if the edge is valid and the tree is traversed further
     */
    virtual bool apply(TreeWalker* walker) const;
};

template<typename Functor>
struct FunctorEdge : public Edge {
  
    FunctorEdge(const Functor& f) : m_functor(f) {};
    
    virtual bool apply(TreeWalker* walker) const {
        
        if(m_functor(walker))
            return Edge::apply(walker);
        
        return false;
    };
    
private:
    Functor m_functor;
};

template<typename Final, typename Constraint, typename Constraint::PimaryOptionType option>
struct ConstraintConditionalEqual {
    
    template<typename Functor>
    struct Type : public Edge {
        
        Type(const Functor& f) : m_functor(f) {};
        
        virtual bool apply(TreeWalker* walker) const {
            
            auto cwalker = static_cast<ConstraintWalker<typename Final::Kernel>*>(walker);
            auto cons = cwalker->template getConstraint<Constraint>(Final::template constraintIndex<Constraint>::value);
            if(!cons)
                return false;
            
            if(cons->template getOption<0>() != option)
                return false;
            
            m_functor(cons->getPrimitiveConstraint());
                            
            return true;
        };
        
    private:
        Functor m_functor;
    };
};
    
struct Node {
  
    virtual ~Node() {    
        //we own the edges,  hence we are responsible for deleting them
        for (Edge* e : m_edges)
            delete e;
    };
    
    /**
     * @brief Connect this node to annother via a specific Edge
     *
     * This function can be used to connect a node to annother one. As the connections are always 
     * conditional the Edge evaluating the condition needs to be supplied.
     * \Note The reduction tree is a static tree, hence there is no identifier for the made connection
     * and also no way to disconnect again
     *
     * \param node The node we build a connection to
     * \param edge The GeometryEdge which evaluates the condition for the transition
     */   
    template<typename Functor>
    void connect(Node* node, Functor func);    
    
    template<typename CEdge, typename Functor> 
    void connectConditional(Node* node, Functor func) {
        typedef typename CEdge::template Type<Functor> EdgeType;
        connect(node, new EdgeType(func));
    };
    
    /**
     * @brief Further traverse the tree
     * 
     * This function is used to execute the node. The base version provided here has no action, hence
     * all it does is to anaylse all edges if they are able to process the walker and apply the walker
     * to a edge if possible. For any actions that should be done by the node a new class must be
     * derived which overrides this method with custom behavior. 
     * \Note When overriding this function you must call it after imüplementing the custom behavior to 
     * ensure the tree traversal
     * 
     * @return bool True if an edge for further traversing was found, false otherwise
     */
    bool apply(TreeWalker* walker) {

        walker->setFinalNode(this);
        for (Edge* e : m_edges) {
            if (e->apply(walker)) {
                return true;
                break;
            }
        };
        return false;
    }; 

private:
    std::vector<Edge*> m_edges;
};

template<typename Functor>
void Node::connect(Node* node, Functor func) {

    auto edge = new FunctorEdge<Functor>(func);
    connect(node, static_cast<Edge*>(edge));
};

template<>
inline void Node::connect<Edge*>(Node* node, Edge* edge) {

    edge->source = this;
    edge->target  = node;
    m_edges.push_back(edge);
};

//we can create the apply function of Edge only now as node must be fullydefined...
inline bool Edge::apply(TreeWalker* walker) const {
    
    target->apply(walker);
    return true;
};

/**
 * @brief Node for Geometries
 * 
 * This node is used for reduction of constraints, it provides an extra functions to build equations
 * for the Geometry it represents. 
 */
template<typename Kernel>
struct GeometryNode : public Node {
    
    typedef shedule::FlowGraph::Node  FlowNode;
    typedef std::shared_ptr<numeric::Calculatable<Kernel>>  Equation;

    virtual Equation buildGeometryEquation(TreeWalker* walker) const = 0;

    virtual std::pair< Equation, FlowNode>
    buildGeometryEquationNode(TreeWalker* walker, shedule::FlowGraph& flowgraph) const = 0;
};

/**
 * @brief Node for Undependend Geometry
 * 
 * This node is used as starting point in the resuction tree, it provides a plain undependend Geometry.
 * As we know the numeric geometry type we have to use we can live only with the primitive Geometry,
 * no need to know any derived class.
 * 
 * \tparam Kernel The mathematical kernel use dto initialize the primitive geometry
 * \tparam G      The primitive geometry to use in the calculation
 */
template<typename Kernel, template<class> class G>
struct UndependendGeometryNode : public GeometryNode<Kernel> {
   
    typedef typename GeometryNode<Kernel>::FlowNode FlowNode;
    typedef typename GeometryNode<Kernel>::Equation Equation;
    
    Equation buildGeometryEquation(TreeWalker* walker) const override {
        TargetWalker<Kernel, G>* gwalker = static_cast<TargetWalker<Kernel, G>*>(walker);
        //create the new primitive geometry and set the initial value
        auto geom = std::make_shared<numeric::Geometry<Kernel, G>>();       
        geom->output() = gwalker->getTargetPrimitive();
        
        return geom;
    };
    
    std::pair< Equation, FlowNode > 
    buildGeometryEquationNode(TreeWalker* walker, shedule::FlowGraph& flowgraph) const override {
        
        TargetWalker<Kernel, G>* gwalker = static_cast<TargetWalker<Kernel, G>*>(walker);
        //create the new primitive geometry and set the initial value
        auto geom = std::make_shared<numeric::Geometry<Kernel, G>>();       
        geom->output() = gwalker->getTargetPrimitive();
        
        return std::make_pair(geom, flowgraph.newActionNode([=](const shedule::FlowGraph::ContinueMessage& m){
            geom->calculate();
        }));
    };
};

/**
 * @brief Node for derived Geometry
 * 
 * This node handles derived geometry in the reduction tree. It takes care that the correct numeric 
 * geometry is created and that all input equations are connected correctly. 
 */
template<typename DerivedG>
struct DependendGeometryNode : public GeometryNode<typename DerivedG::KernelType> {

    typedef typename DerivedG::KernelType   Kernel;
    typedef typename DerivedG::OutputType   Output;
    typedef typename DerivedG::InputType    Input;
    
    typedef typename GeometryNode<Kernel>::FlowNode FlowNode;
    typedef typename GeometryNode<Kernel>::Equation Equation;
    
    virtual Equation buildGeometryEquation(TreeWalker* walker) const {
        auto* gwalker = static_cast<TargetWalker<Kernel, geometry::extractor<Output>::template primitive>*>(walker);
        //create the new primitive geometry and set the initial value
        auto geom = std::make_shared<DerivedG>();       
        geom->output() = gwalker->getTargetPrimitive();
        
        //set the input for our unary equation
        auto eqn = std::static_pointer_cast<numeric::Equation<Kernel, Input>>(gwalker->getInputEquation());
        geom->setInputEquation(eqn);
        geom->takeInputOwnership(true);
        
        return geom;
    }
    
    virtual std::pair< Equation, FlowNode > buildGeometryEquationNode(TreeWalker* walker, shedule::FlowGraph& flowgraph) const {
        
        auto* gwalker = static_cast<TargetWalker<Kernel, geometry::extractor<Output>::template primitive>*>(walker);
        //create the new primitive geometry and set the initial value
        auto geom = std::make_shared<DerivedG>();       
        geom->output() = gwalker->getTargetPrimitive();
        
        //set the input for our unary equation
        auto eqn = std::static_pointer_cast<numeric::Equation<Kernel, Input>>(gwalker->getInputEquation());
        geom->setInputEquation(eqn);
        geom->takeInputOwnership(true);
        
        return std::make_pair(geom, flowgraph.newActionNode([=](const shedule::FlowGraph::ContinueMessage& m){
            geom->calculate();
        }));
    };
};

struct EdgeReductionTree {

    /**
     * @brief Analyses the global edges and finds the best reduction result
     *
     * The caller owns the returned TreeWalker pointer.
     * 
     * @remark The function is reentrant but is not safe to be called on the same data from multiple threads
     *
     * @param g The graph the local edge belongs to
     * @param e The local edge to analyse for reduction
     * @return void
     */
    virtual reduction::TreeWalker* apply(symbolic::Geometry* source, symbolic::Geometry* target,
                                         std::vector<symbolic::Constraint*> constraints) = 0;
};

template<typename Kernel, template<class> class SourceGeometry, template<class> class TargetGeometry>
struct GeometryEdgeReductionTree : public EdgeReductionTree {

    typedef typename Kernel::Scalar Scalar;

    GeometryEdgeReductionTree() {
        //create the tree source note, it must always be available
        m_sourceNode = new reduction::UndependendGeometryNode<Kernel, TargetGeometry>();
    };
    
    virtual ~GeometryEdgeReductionTree() {
        //we own all nodes, delete them!      
        delete m_sourceNode;
        for(auto& iter : m_nodesMap) {
            if(iter.second) 
                delete iter.second;
            
        }
    };
    
    reduction::Node* getSourceNode() {return m_sourceNode;};
    
    template<typename T>
    reduction::Node* getTreeNode() { 
        auto iter = m_nodesMap.find(std::type_index(typeid(T)));
        if(iter != m_nodesMap.end())
            return iter->second;
        
        m_nodesMap[std::type_index(typeid(T))] = new DependendGeometryNode<T>();
        return m_nodesMap[std::type_index(typeid(T))];
    };


    virtual reduction::TreeWalker* apply(symbolic::Geometry* source, symbolic::Geometry* target,
                                         std::vector<symbolic::Constraint*> constraints) {

        dcm_assert(source != target);
        //dcm_assert(dynamic_cast<TypeGeometry<Kernel, SourceGeometry>*>(source) != NULL);
        //dcm_assert(dynamic_cast<TypeGeometry<Kernel, SourceGeometry>*>(target) != NULL);

        //get the primitive geometries
        const SourceGeometry<Kernel>& pg1 = static_cast<symbolic::TypeGeometry<Kernel, SourceGeometry>*>(source)->getPrimitveGeometry();
        const TargetGeometry<Kernel>& pg2 = static_cast<symbolic::TypeGeometry<Kernel, TargetGeometry>*>(target)->getPrimitveGeometry();

        //create a new treewalker and set it up
        auto walker = new SourceTargetWalker<Kernel, SourceGeometry, TargetGeometry>(pg1, pg2);
        walker->setInitialNode(getSourceNode());
        walker->setConstraintPool(constraints);
        
        //start the calculation
        getSourceNode()->apply(walker);
        return walker;
    };
    
protected:
    reduction::UndependendGeometryNode<Kernel, TargetGeometry>* m_sourceNode;
    std::unordered_map<std::type_index, reduction::Node*>       m_nodesMap;
};

} //reduction

namespace symbolic {

template<typename Kernel>
struct Reduction : public numeric::EquationBuilder<Kernel> {
    
    //the general reduction result: how many parameters can be eliminated
    int parameterReduction() {
        return m_paramReduction;
    }

private:
    int m_paramReduction; // the amount of degrees of freedom reduced by this reduction
};

template<typename Kernel>
struct ConstraintReduction : public Reduction<Kernel> {
    
    typedef shedule::FlowGraph::Node  Node;
    typedef std::shared_ptr<numeric::Calculatable<Kernel>> CalcPtr;
    
    ConstraintReduction(const graph::LocalVertex& source, const graph::LocalVertex& target,
                        reduction::ConstraintWalker<Kernel>* sw, reduction::ConstraintWalker<Kernel>* tw) 
        : m_sourceWalker(sw), m_targetWalker(tw), m_sourceVertex(source), m_targetVertex(target) {};
    
    virtual ~ConstraintReduction() {
        //we own the walkers
        delete m_sourceWalker;
        delete m_targetWalker;
    }
        
    //this function is used to create a default numeric geometry for the given symbolic one.
    virtual CalcPtr createGeometry(const graph::LocalVertex& vertex) {
               
        if(vertex == m_sourceVertex) {
            auto node = static_cast<reduction::GeometryNode<Kernel>*>(m_sourceWalker->getInitialNode());
            return node->buildGeometryEquation(m_sourceWalker);
        }
        else  {
            auto node = static_cast<reduction::GeometryNode<Kernel>*>(m_targetWalker->getInitialNode());
            return node->buildGeometryEquation(m_targetWalker);
        }
    };

    virtual std::pair< CalcPtr, Node>
    createGeometryNode(const graph::LocalVertex& vertex, shedule::FlowGraph& flowgraph) {
        
        if(vertex == m_sourceVertex) {
            auto node = static_cast<reduction::GeometryNode<Kernel>*>(m_sourceWalker->getInitialNode());
            return node->buildGeometryEquationNode(m_sourceWalker, flowgraph);
        }
        else  {
            auto node = static_cast<reduction::GeometryNode<Kernel>*>(m_targetWalker->getInitialNode());
            return node->buildGeometryEquationNode(m_targetWalker, flowgraph);
        }
    };
    
    //this function creates a reduced numeric geometry equation. For this the equation it depends on 
    //is required
    virtual CalcPtr createReducedGeometry(const graph::LocalVertex& vertex, CalcPtr basegeometry) {
        
        if(vertex == m_sourceVertex) {
            auto node = static_cast<reduction::GeometryNode<Kernel>*>(m_sourceWalker->getFinalNode());
            return node->buildGeometryEquation(m_sourceWalker);
        }
        else  {
            auto node = static_cast<reduction::GeometryNode<Kernel>*>(m_targetWalker->getFinalNode());
            return node->buildGeometryEquation(m_targetWalker);
        }
    };
     
    //this function is used to create a reduced numeric geometry equation and a node for it in the 
    //corresponding flow graph. For this the equation it depends on is required
    virtual std::pair< CalcPtr, Node>
    createReducedGeometryNode(const graph::LocalVertex& vertex, CalcPtr basegeometry,
                              shedule::FlowGraph& flowgraph) {
        
        if(vertex == m_sourceVertex) {
            auto node = static_cast<reduction::GeometryNode<Kernel>*>(m_sourceWalker->getFinalNode());
            return node->buildGeometryEquationNode(m_sourceWalker, flowgraph);
        }
        else  {
            auto node = static_cast<reduction::GeometryNode<Kernel>*>(m_targetWalker->getFinalNode());
            return node->buildGeometryEquationNode(m_targetWalker, flowgraph);
        }
    };
    
    //this function is used to create all default constraint equations
    virtual std::vector<CalcPtr>
    createEquations(CalcPtr g1, CalcPtr g2) {
        
    };
                       
    //this function is used to create a node in the calculation flow graph with all remaining
    //constraint Equations
    virtual std::pair<std::vector<CalcPtr>, Node>
    createEquationsNode(CalcPtr g1, CalcPtr g2, shedule::FlowGraph& flow) {
    
        
    };

    //this function is used to create all default constraint equations
    virtual std::vector<CalcPtr>
    createReducedEquations(CalcPtr g1, CalcPtr g2) {
        
    };
                       
    //this function is used to create a node in the calculation flow graph with all remaining
    //constraint Equations
    virtual std::pair<std::vector<CalcPtr>, Node>
    createReducedEquationsNode(CalcPtr g1, CalcPtr g2, shedule::FlowGraph& flow) {
        
    };
    
private:
    graph::LocalVertex                   m_sourceVertex, m_targetVertex;
    reduction::ConstraintWalker<Kernel>* m_sourceWalker; //reduction result for the source geometry
    reduction::ConstraintWalker<Kernel>* m_targetWalker; //reduction result for the target geometry 
    std::vector<symbolic::Constraint*>   m_constraints;  //all constraints in the edge
};
    
template<typename Final>
struct Reducer {
    
    typedef typename Final::Kernel Kernel;
    
    Reducer() {
    
        int size = mpl::size<typename Final::GeometryList>::type::value;
        m_treeArray.resize(boost::extents[size][size]);
        
        //build up the default reduction nodes
        //mpl trickery to get a sequence counting from 0 to the size of stroage entries
        typedef mpl::range_c<int,0,
                mpl::size<typename Final::GeometryList>::value> StorageRange;
        //now iterate that sequence so we can access all storage elements with knowing the position
        //we are at
        utilities::RecursiveSequenceApplyer<typename Final::GeometryList, ReductionTreeCreator> r(m_treeArray);
        mpl::for_each<StorageRange>(r);
    };
    
    void reduce(std::shared_ptr<typename Final::Graph> g, graph::LocalEdge edge) {
        
        //get the geometry used in this edge
        symbolic::Geometry* source = g->template getProperty<symbolic::GeometryProperty>(g->source(edge));
        symbolic::Geometry* target = g->template getProperty<symbolic::GeometryProperty>(g->target(edge));
        
        //get the two reduction trees for this geometry combination
        reduction::EdgeReductionTree* stTree = m_treeArray[source->type][target->type];
        reduction::EdgeReductionTree* tsTree = m_treeArray[target->type][source->type];
        
        //get all constraints
        std::vector<symbolic::Constraint*> constraints;
        typedef typename Final::Graph::global_edge_iterator iterator;
        std::pair<iterator, iterator> it = g->getGlobalEdges(edge);
        for (; it.first != it.second; ++it.first)
            constraints.push_back(g->template getProperty<symbolic::ConstraintProperty>(*it.first));

        //calculate both results
        auto* stWalker = static_cast<reduction::ConstraintWalker<Kernel>*>(stTree->apply(source, target, constraints));
        auto* tsWalker = static_cast<reduction::ConstraintWalker<Kernel>*>(stTree->apply(target, source, constraints));
        
        //build the reduction and set store it in the graph. Make sure any pointer already stored is
        //deleted properly, especially the walkers
        auto reduction = g->template getProperty<numeric::EquationBuilderProperty<Kernel>>(edge);
        if(reduction)
            delete reduction;
        
        reduction = new ConstraintReduction<typename Final::Kernel>(g->source(edge), g->target(edge), tsWalker, stWalker);
        g->template setProperty<numeric::EquationBuilderProperty<Kernel>>(edge, reduction);
    };
    
private:
    boost::multi_array<reduction::EdgeReductionTree*,2> m_treeArray;
       
    template<typename Sequence>
    struct ReductionTreeCreator {
    
        boost::multi_array<reduction::EdgeReductionTree*,2>& m_treeArray;
        
        ReductionTreeCreator(boost::multi_array<reduction::EdgeReductionTree*,2>& r) 
            : m_treeArray(r) {};
            
        template<typename N1, typename N2>
        void operator()() {
        
            typedef typename mpl::at<Sequence, N1>::type t1;
            typedef typename mpl::at<Sequence, N2>::type t2;
            
            auto node1 = new reduction::GeometryEdgeReductionTree<Kernel, 
                                geometry::extractor<t1>::template primitive,
                                geometry::extractor<t2>::template primitive >();
                                
            auto node2 = new reduction::GeometryEdgeReductionTree<Kernel, 
                                geometry::extractor<t2>::template primitive,
                                geometry::extractor<t1>::template primitive >();
        
            int idx1 = Final::template initGeometryIndex<t1>::value;
            int idx2 = Final::template initGeometryIndex<t2>::value;
            
            m_treeArray[idx1][idx2] = node1;
            m_treeArray[idx2][idx1] = node2;
        };
    };
};

}//symbolic
}//dcm

#endif //DCM_ANALYSE_H

