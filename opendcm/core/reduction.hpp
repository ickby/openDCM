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
namespace symbolic {
   

/**
 * @brief Structure to setup the numeric solver from reduction results
 *
 * This struct is the entry point for the numeric solver setup. It has all the functions needed to get the
 * numeric geometries and constraints.
 */
template<typename Kernel>
struct Reduction {
    
    typedef shedule::FlowGraph::Node  Node;
    typedef std::shared_ptr<numeric::Calculatable<Kernel>> CalcPtr;
    
    //returns how many parameters are saved when the reduced version is used
    int parameterReduction() {
        return m_paramReduction;
    }
    
    //this function is used to create a default numeric geometry for the given symbolic one.
    virtual CalcPtr createGeometry(symbolic::Geometry& geometry) = 0;
     
    //this function creates a default numeric geometry equation and a node for it in the 
    //corresponding flow graph
    virtual std::pair< CalcPtr, Node>
    createGeometryNode(symbolic::Geometry& geometry, shedule::FlowGraph& flowgraph) = 0;
    
    //this function creates a reduced numeric geometry equation. For this the equation it depends on 
    //is required
    virtual CalcPtr createReducedGeometry(symbolic::Geometry& geometry, CalcPtr basegeometry) = 0;
     
    //this function is used to create a reduced numeric geometry equation and a node for it in the 
    //corresponding flow graph. For this the equation it depends on is required
    virtual std::pair< CalcPtr, Node>
    createReducedGeometryNode(symbolic::Geometry& geometry, CalcPtr basegeometry,
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
    
protected:
    CalcPtr              m_geometryEquation;
    std::vector<CalcPtr> m_residualEquations;
    
    int m_paramReduction; // the amount of degrees of freedom reduced by this reduction
};

template<typename Kernel>
struct ResultProperty {
    typedef Reduction<Kernel>* type;
    struct default_value {
        Reduction<Kernel>* operator()() {
            return nullptr;
        };
    };
};

/**
 * @brief The tree structure used for Constrain reduction
 * 
 * 
 */
namespace reduction {

//TODO: Currently the geometry equations are created during traversal. This is highly inefficient: 
//      for one it means creating and destroying many unused equations during traversal, and on the
//      other side means when creating a flowgraph node only the virtual execute is available, 
//      resulting in 2 virtual calls for one execution instead of minimal 1. Solution is simple: 
//      split node action from traversal and store the node in the TreeWalker. Than the walker can 
//      call the creation method of the node when requested and with a flowgraph argument
    
struct Edge;

/**
 * @brief Data structure for tree traversal
 * 
 * As the reduction tree should be traversed in parallel the nodes and edges must be reentrant, hence
 * they do not hold any data and only provide algorithms. The data for the traversal is stored in the
 * treewalker.
 */
struct TreeWalker {
    
};

struct Node;

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

template<typename Constraint, typename Constraint::
    
struct Node {
  
    ~Node() {      
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
    virtual bool apply(TreeWalker* walker) {
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
void Node::connect<Edge*>(Node* node, Edge* edge) {

    edge->source = this;
    edge->target  = node;
    m_edges.push_back(edge);
};

//we can create the apply function of Edge only now as node must be fullydefined...
bool Edge::apply(TreeWalker* walker) const {
        
    target->apply(walker);
    return true;
};

template<typename Kernel, template<class> class Primitive>
struct GeometryWalker : public TreeWalker {
    
    typedef std::shared_ptr<numeric::Equation<Kernel, Primitive<Kernel>>> Geometry;
    typedef std::shared_ptr<numeric::Calculatable<Kernel>>                Equation;
    
    GeometryWalker(const Primitive<Kernel>& prim) : m_primitive(prim) {};
 
    const Primitive<Kernel>& getPrimitive() {return m_primitive;};
    
    void        setGeometry(Geometry g) {m_geometry = g;};
    Geometry    getGeometry() {return m_geometry;};
    
    void        setCummulativeInputEquation(Equation e) {m_inputEqn = e;};
    Equation    getCummulativeInputEquation() {return m_inputEqn;};
    
private:
    const Primitive<Kernel>&    m_primitive; //primitive holding the value
    Equation                    m_inputEqn;  //cummulative input equation
    Geometry                    m_geometry;  //numeric geometry 
};

template<typename Kernel, template<class> class SourcePrimitive, template<class> class TargetPrimitive>
struct ConstraintWalker : public GeometryWalker<Kernel, TargetPrimitive> {
       
    typedef typename std::vector<symbolic::Constraint*>::const_iterator iterator;
    
    ConstraintWalker(const SourcePrimitive<Kernel>& sprim, const TargetPrimitive<Kernel>& tprim) :
                 GeometryWalker<Kernel, TargetPrimitive>(tprim), m_sourcePrimitive(sprim) {};
 
    const SourcePrimitive<Kernel>& getSourcePrimitive() {return m_sourcePrimitive;};
    
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
    const SourcePrimitive<Kernel>&      m_sourcePrimitive; //primitive holding the value of the source
    std::vector<symbolic::Constraint*>  m_constraintPool;  //all the constraints we want to reduce
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
struct GeometryNode : public Node {

    virtual bool apply(TreeWalker* walker) {
    
        GeometryWalker<Kernel, G>* gwalker = static_cast<GeometryWalker<Kernel, G>*>(walker);
        
        //create the new primitive geometry and set the initial value
        auto geom = std::make_shared<numeric::Geometry<Kernel, G>>();       
        geom->output() = gwalker->getPrimitive();
        
        //set the new value in the walker for further processing
        gwalker->setGeometry(geom);
        
        //further traverse the tree
        Node::apply(walker);

        return true;
    }
};

/**
 * @brief Node for derived Geometry
 * 
 * This node handles derived geometry in the reduction tree. It takes care that the correct numeric 
 * geometry is created and that all input equations are connected correctly. 
 */
template<typename DerivedG>
struct DependendGeometryNode : public Node {

    virtual bool apply(TreeWalker* walker) {
        
        typedef typename DerivedG::KernelType   Kernel;
        typedef typename DerivedG::OutputType   Output;
        typedef typename DerivedG::InputType    Input;
                
        auto* gwalker = static_cast<GeometryWalker<Kernel, geometry::extractor<Output>::template primitive>*>(walker);
        
        //create the new primitive geometry and set the initial value
        auto geom = std::make_shared<DerivedG>();       
        geom->output() = gwalker->getPrimitive();
        
        //set the input for our unary equation
        auto eqn = std::static_pointer_cast<numeric::Equation<Kernel, Input>>(gwalker->getCummulativeInputEquation());
        geom->setInputEquation(eqn);
        geom->takeInputOwnership(true);
        
        //set the new value in the walker for further processing
        //note that the walker has the only shared_ptr of the equation, hence if it is a unary 
        //equation which has ownership over the input it will be deleted and hence no double 
        //ownership will occure
        gwalker->setGeometry(geom);
        
        //further traverse the tree
        Node::apply(walker);
        
        return true;
    }
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
        m_sourceNode = new reduction::GeometryNode<Kernel, TargetGeometry>();
    };
    
    virtual ~GeometryEdgeReductionTree() {
        //we own all nodes, delete them!        
        for(auto iter = m_nodesMap.begin(); iter != m_nodesMap.end(); ++iter)
            delete iter->second;
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
        const SourceGeometry<Kernel>& pg1 = static_cast<TypeGeometry<Kernel, SourceGeometry>*>(source)->getPrimitveGeometry();
        const TargetGeometry<Kernel>& pg2 = static_cast<TypeGeometry<Kernel, SourceGeometry>*>(target)->getPrimitveGeometry();

        //create a new treewalker and set it up
        auto walker = new ConstraintWalker<Kernel, SourceGeometry, TargetGeometry>(pg1, pg2);
        walker->setConstraintPool(constraints);
        
        //start the calculation
        getSourceNode()->apply(walker);
        return walker;
    };
    
protected:
    reduction::GeometryNode<Kernel, TargetGeometry>*      m_sourceNode;
    std::unordered_map<std::type_index, reduction::Node*> m_nodesMap;
};

} //reduction

template<typename Final>
struct Reducer {
    
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
    
    void reduce(typename Final::Graph& g, graph::LocalEdge edge) {
        
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
        reduction::TreeWalker* stWalker = stTree->apply(source, target, constraints);
        reduction::TreeWalker* tsWalker = stTree->apply(target, source, constraints);
        
        //build the reduction and set store it in the graph. Make sure any pointer already stored is
        //deleted properly, especially the walkers
    };
    
private:
    boost::multi_array<symbolic::reduction::EdgeReductionTree*,2> m_treeArray;
       
    template<typename Sequence>
    struct ReductionTreeCreator {
    
        boost::multi_array<symbolic::reduction::EdgeReductionTree*,2>& m_treeArray;
        
        ReductionTreeCreator(boost::multi_array<symbolic::reduction::EdgeReductionTree*,2>& r) 
            : m_treeArray(r) {};
            
        template<typename N1, typename N2>
        void operator()() {
        
            typedef typename mpl::at<Sequence, N1>::type t1;
            typedef typename mpl::at<Sequence, N2>::type t2;
            
            auto node1 = new symbolic::reduction::GeometryEdgeReductionTree<Final, 
                                geometry::extractor<t1>::template primitive,
                                geometry::extractor<t2>::template primitive >();
                                
            auto node2 = new symbolic::reduction::GeometryEdgeReductionTree<Final, 
                                geometry::extractor<t2>::template primitive,
                                geometry::extractor<t1>::template primitive >();
        
            int idx1 = Final::template geometryIndex<t1>::value;
            int idx2 = Final::template geometryIndex<t2>::value;
            
            m_treeArray[idx1][idx2] = node1;
            m_treeArray[idx2][idx1] = node2;
        };
    };
};

}//symbolic
}//dcm

#endif //DCM_ANALYSE_H

