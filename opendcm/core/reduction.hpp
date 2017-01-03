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

#ifndef DCM_REDUCTION_H
#define DCM_REDUCTION_H

#include <Eigen/Core>

#include "geometry.hpp"
#include "clustergraph.hpp"
#include "constraint.hpp"
#include "scheduler.hpp"
#include "utilities.hpp"

#include <boost/tuple/tuple.hpp>
#include <boost/iterator/zip_iterator.hpp>
#include <boost/multi_array.hpp>
#include <unordered_map>
#include <typeindex>
#include <vector>
                

namespace dcm {
namespace numeric {   

/**
 * @brief Structure to setup the numeric equations for a solver at a graph edge
 *
 * This struct is the entry point for the numeric solver setup. It has all the functions needed to get the
 * numeric geometries and constraints for a \ref ClusterGraph edge. As the usage of vertices as un- or 
 * dependend geometry is not yet decided at creation time this clas offers the possibility to create both, 
 * dependend on the callers needs. The caller is responsible for calling the appropriate functions for 
 * reduced/unreduced geometry and constraints and to not mix them. Futhermore the class provides functions
 * to create equations directly in a \ref FlowGraph which allows for more efficient calculation than adding 
 * the equation to a \ref FlowGraph after the creation. When using the provided methods a minimal amount of
 * virtual calls is assured.
 * 
 * @note The returned equations are not yet initialised and hence cannot be used directly. It is the
 * callers responsibility to call \ref init for all equations
 */
template<typename Kernel>
struct EquationBuilder {
    
    virtual ~EquationBuilder(){};
    
    typedef shedule::FlowGraph::Node  Node;
    typedef std::shared_ptr<numeric::Calculatable<Kernel>> CalcPtr;
    
    
    /**
     * @brief Create a undependend geometry equation 
     * 
     * This function creates the undependend geometry equation for a given vertex. This vertex must 
     * be the source or target of the edge where the EquationBuilder is added to. The result does 
     * ignore all reduction results done on the edge.
     * 
     * @param vertex The vertex that describes the geometry in the clustergraph
     * @return CalcPtr The numeric eqution
     */
    virtual CalcPtr createGeometry(const graph::LocalVertex& vertex) = 0;
     
    /**
     * @brief Create a undependend geometry equation in a \ref FlowGraph
     * 
     * This function creates the undependend geometry equation for a given vertex directly in a 
     * \ref FlowGraph The given vertex must be the source or target of the edge where the 
     * EquationBuilder is added to. The result does ignore all reduction results done on the edge. 
     * The returned node is unconnected. To allow initialisation of the equation it creates it is 
     * provided too.
     * 
     * @param vertex The vertex that describes the geometry in the clustergraph
     * @param flowgraph The \ref FlowGraph the node should be created in
     * @return td::pair< CalcPtr, Node> The numeric eqution and the Node in the \ref FlowGraph
     */
    virtual std::pair< CalcPtr, Node>
    createGeometryNode(const graph::LocalVertex& vertex, shedule::FlowGraph& flowgraph) = 0;
    
    /**
     * @brief Create a dependend geometry equation
     * 
     * This function creates the dependend geometry equation for a given vertex. This vertex must 
     * be the source or target of the edge where the EquationBuilder is added to. The created equation 
     * is a \ref InputEquation depending on the edges other vertex geometry, which is not provided, and
     * the constraints which have been used for the reduction. Hence if this EquationBuilder is
     * added to a certain edge, this function called with the target vertex creates a geometry of the
     * target type dependend on the source type and the constraints of this edge.
     * @param vertex The vertex of the depending geometry to create
     * @return CalcPtr The numeric eqution
     */
    virtual CalcPtr createReducedGeometry(const graph::LocalVertex& vertex) = 0;
     
    
    /**
     * @brief Creates a dependend geometry equation in a \ref FlowGraph
     * 
     * This function does the same as \ref createReducedGeometry but adds the equation directly in a
     * \ref FlowGraph and ensures the most efficient calculation of the node by reducing virtual 
     * function calls. The returned node is unconnected. To allow initialisation of the equation it
     * creates it is provided too.
     * @param vertex The vertex of the depending geometry to create
     * @param flowgraph The \ref FlowGraph the node should be created in
     * @return std::pair< CalcPtr, Node > The numeric eqution and the Node in the \ref FlowGraph
     */
    virtual std::pair< CalcPtr, Node>
    createReducedGeometryNode(const graph::LocalVertex& vertex, shedule::FlowGraph& flowgraph) = 0;
    

    /**
     * @brief Creates numeric equations for all available constraints
     * 
     * This function creates numeric equations for all available constraints, no matter what the 
     * reduction did result in. The created equations of type BinaryEquation have the geometries 
     * as input, which are at the EquationBuilders edge source/target. The geometries at source and target vertex
     * of the edge should be used as source and target equation for this function, mixing those will
     * lead to faulty equation claculation! 
     * 
     * @param sourceGeometry Geometry equation at the edges source vertex
     * @param targetGeometry Geometry equation at the edges target vertex
     * @return std::vector< CalcPtr > The vector with BinaryEquation of all constraints 
     */    
    virtual std::vector<CalcPtr>
    createEquations(CalcPtr sourceGeometry, CalcPtr targetGeometry) = 0;
                       

    /**
     * @brief Creates numeric equations for all available constraints directly in a \ref FlowGraph
     * 
     * This method does exactly the same as \ref createEquations but creates those equations in the 
     * provided \ref FlowGraph. It is ensured that the calculation of the created equtions in the 
     * graph node is as efficient as possible. The returned node is unconnected. To allow initialisation 
     * of the equations created they are provided too.
     * 
     * @param sourceGeometry Geometry equation at the edges source vertex
     * @param targetGeometry Geometry equation at the edges target vertex
     * @param flow The \ref FlowGraph the node should be created in
     * @return std::pair< std::vector< CalcPtr >, Node > The constraint BinaryEquation vector and the Node 
     *                                                   in the \ref FlowGraph
     */
    virtual std::pair<std::vector<CalcPtr>, Node>
    createEquationsNode(CalcPtr sourceGeometry, CalcPtr targetGeometry, shedule::FlowGraph& flow) = 0;

    /**
     * @brief Creates numeric equations the remaining constraints in a \ref FlowGraph
     * 
     * This function creates numeric equations for the remaining constraints after the reduction 
     * process. The created equations of type BinaryEquation have the geometries 
     * as input, which are at the EquationBuilders edge source/target. The geometries at source and 
     * target vertex of the edge should be used as source and target equation for this function, mixing 
     * those will lead to faulty equation claculation! 
     * \note It is possible that there are no remaining equations and the returned vector is empty
     * 
     * @param sourceGeometry Geometry equation at the edges source vertex
     * @param targetGeometry Geometry equation at the edges target vertex
     * @return std::vector< CalcPtr > The vector with BinaryEquation of all remaining constraints 
     */  
    virtual std::vector<CalcPtr>
    createReducedEquations(const graph::LocalVertex& target, CalcPtr sourceGeometry, CalcPtr targetGeometry) = 0;
                       
    /**
     * @brief Creates numeric equations for all available constraints directly in a \ref FlowGraph
     * 
     * This method does exactly the same as \ref createReducedEquations but creates those equations in the 
     * provided \ref FlowGraph. It is ensured that the calculation of the created equtions in the 
     * graph node is as efficient as possible. The returned node is unconnected. To allow initialisation 
     * of the equations created they are provided too.
     * \note It is possible that there are no remaining equations and the returned vector is empty. In
     * this case the created Node should be deleted again
     * 
     * @param sourceGeometry Geometry equation at the edges source vertex
     * @param targetGeometry Geometry equation at the edges target vertex
     * @param flow The \ref FlowGraph the node should be created in
     * @return std::pair< std::vector< CalcPtr >, Node > The constraint BinaryEquation vector and the Node 
     *                                                   in the \ref FlowGraph
     */
    virtual std::pair<std::vector<CalcPtr>, Node>
    createReducedEquationsNode(const graph::LocalVertex& target, CalcPtr sourceGeometry, 
                               CalcPtr targetGeometry, shedule::FlowGraph& flow) = 0;
};

/**
 * @brief Store a EquationBuilder at a graph edge
 *
 * EquationBuilder class is used to convert a symbolic edge with its symbolic vertices to numeric 
 * equations. This builder is the result of an symbolic analysis and must be stored at the edge which
 * was analysed. This property does make this possible, it is added as a LocalEdge property to the 
 * ClusterGraph. It stores the EquationBuilder as pointer to allow polymoprphism. The initial value 
 * is a nullptr, hence always check for validity! 
 * 
 * \tparam Kernel The Kernel used 
 **/
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
 * @brief The tree structure used for constrain reduction
 * 
 * 
 */
namespace reduction {
   
struct Connection;
struct Node;

/**
 * @brief Data structure for tree traversal
 * 
 * As the reduction tree should be traversed in parallel the nodes and edges must be reentrant, hence
 * they must not hold any data and only provide algorithms. The data for the traversal is stored in the
 * treewalker. This base class has the sole purpose of recording at which node the traversal started 
 * and at which it ended. For everything else, data storage and processing, derived classes are 
 * responsible
 */
struct TreeWalker {
    
    /** @brief Retrieve node where the tree traversal started
     * @return dcm::reduction::Node* Pointer to initial node
     */
    Node* getInitialNode()      {return m_initialNode;};
    
    /** @brief Set initial node before starting the traversal
     */    
    void  setInitialNode(Node* n) {m_initialNode = n;};

    /**
     * @brief Retrieve the node at which the tree traversal ended     * 
     * @return dcm::reduction::Node* Pointer to final node
     */
    Node* getFinalNode()      {return m_finalNode;};
    
    /**
     * @brief Set the node at which the walker currently is during traversal
     */
    void  setCurrentNode(Node* n) {m_finalNode = n;};
    
private:
    Node* m_initialNode;
    Node* m_finalNode;
};

/**
 * @brief Tree walker that stores constraints which can be used for reduction
 * 
 * This walker stores all constraints which can be used as edge indicators to traverse the tree. It
 * allows the edges to access the available constraints for decision if they are a valid connection.
 * If an endge uses a constraint it can remove it from the pool of available constraints.
 * The class exposes a container interface for the symbolic constraints and hence can be used in 
 * range based for loops.
 * 
 * @tparam Kernel The math kernel in use
 */
template<typename Kernel>
struct ConstraintWalker : public TreeWalker {
    
    typedef std::shared_ptr<numeric::Calculatable<Kernel>>                            Equation;
    typedef std::tuple<symbolic::Constraint*,symbolic::Geometry*,symbolic::Geometry*> SymbolicTuple;
    typedef std::vector<SymbolicTuple>                                                SymbolicVector;    
    typedef typename SymbolicVector::const_iterator                                   Iterator;  
    
    /** @brief Fill constraint pool with available symbolic constraints
     */
    void setConstraintPool(const SymbolicVector& vec) {m_constraintPool = vec;};
    
    //interface for accessing constraints
    /** @brief Check if constraint pools is empty      
     * @return bool true if no constraints are available
     */
    bool     empty() const {return m_constraintPool.empty();}
    
    /** @brief Retrieve the amount of available constraints     
     * @return int number of available constraints
     */
    int      size()  const {return m_constraintPool.size();}
    
    /** @brief Iterator to start constraint traversal
     * @return Iterator Iterator which points at the first element in the sequence 
     */
    Iterator begin() const {return m_constraintPool.begin();};
    
    /** @brief Iterator to test end of sequence
     * @return Iterator Iterator which points after the last constraint element
     */
    Iterator end()   const {return m_constraintPool.end();};
    
    /** @brief First constraint in the sequence
     * @return const symbolic::Constraint* First symbolic constraint
     */
    const SymbolicTuple front() const {return m_constraintPool.front();};
    
    /** @brief Last constraint in the sequence
     * @return const symbolic::Constraint* Last symbolic constraint
     */
    const SymbolicTuple back() const {return m_constraintPool.back();};
    
    
    /**
     * @brief Get constraint by type if available
     * 
     * Returns a TypeConstraint if it is available in the sequence. Note that this is a constraint 
     * walker, hence it is only used for reducing edges. It is therefore not allowed to have a 
     * constraint type multiple times in the constraint pool, they would always be redundant or 
     * conflicting. Hence this method is safe to use and there is no danger of not finding a constraint
     * which comes later in the sequence.
     * @tparam T primitive constraint type to retrieve
     * @return symbolic::TypeConstraint<T>* Symbolic TypeConstraint holding the primitive contraint or 
     *                                      nullptr if not available
     */
    template<typename T>
    symbolic::TypeConstraint<T>* getConstraint(int type) {
        for(auto tuple : m_constraintPool) {
            if(std::get<0>(tuple)->type == type)
                return static_cast<symbolic::TypeConstraint<T>*>(std::get<0>(tuple));
        }
        return nullptr;
    };
    
    /**
     * @brief Get constraint tuple including geometries by constraint type
     * 
     * @return const SymbolicTuple& Tuple of constraint and zource/zarget geometry
     */
    SymbolicTuple getConstraintTuple(int ConstraintType) {
    
        for(auto tuple : m_constraintPool) {
            if(std::get<0>(tuple)->type == ConstraintType)
                return tuple;
        }
        return SymbolicTuple();
    };
    
    /**
     * @brief Accept constraint during tree traversal
     * 
     * When a constraint is used by an edge it can be consumed and than not be available anymore in
     * the constraint pool. This method is used for this: when called the given constraint is removed
     * from the pool
     */   
    void acceptConstraint(symbolic::Constraint* c) {
        for(auto tuple : m_constraintPool) {
            if(std::get<0>(tuple) == c) {
                acceptConstraintTuple(tuple);
                break;
            }
        }
    };
    
    /**
     * @brief Accept constraint tuple during tree traversal
     * 
     * When a constraint is used by an edge it can be consumed and than not be available anymore in
     * the constraint pool. This method is used for this: when called the given constraint tuple is removed
     * from the pool
     */ 
    void acceptConstraintTuple(const SymbolicTuple& tuple) {
        m_constraintPool.erase(std::remove(m_constraintPool.begin(), m_constraintPool.end(), tuple), m_constraintPool.end());
    }
    
private:
    SymbolicVector  m_constraintPool;  //all the constraints we want to reduce
};

/**
 * @brief TreeWalker with information about the edges target geometry
 * 
 * This class extends the ConstraintWalker with additional information about the edges target, namely
 * the geometry it holds. It allows to access PrimitiveGeometry stored at the target vertex. As for 
 * this walker it is known that the target is holds a geometry an that the goal is to reduce the 
 * constraints to build a dependend geometrie it also provides an interface to build up the input 
 * equation for this DependendGeometry. It allows to set and access the Equation used as input and 
 * hence allows to build up an stacked equation step by step.
 * 
 * @tparam Kernel The math Kernel in use
 * @tparam Primitive The primitive geometry type stored at the edges source vertex
 */
template<typename Kernel, template<class> class Primitive>
struct TargetWalker : public ConstraintWalker<Kernel> {
    
    typedef ConstraintWalker<Kernel>                                      Inherited;
    typedef std::shared_ptr<numeric::Equation<Kernel, Primitive<Kernel>>> Geometry;
    typedef std::shared_ptr<numeric::Calculatable<Kernel>>                Equation;
    
    TargetWalker(const Primitive<Kernel>& prim) : m_targetPrimitive(prim) {};
 
    /**
     * @brief Access the primitive geometry stored at the edges target vertex
     * @return const Primitive< Kernel >& The stored primitive geometry
     */
    const Primitive<Kernel>& getTargetPrimitive() {return m_targetPrimitive;};
        
    /**
     * @brief Set the commulative input equation
     */
    void        setInputEquation(Equation e) {m_inputEqn = e;};
    
    /**
     * @brief Retrieve the earlyer set input equation
     * @return Equation Can be empty pointer if no equation was set earlyer
     */    
    Equation    getInputEquation() {return m_inputEqn;};
    
private:
    const Primitive<Kernel>&    m_targetPrimitive; //primitive holding the value
    Equation                    m_inputEqn;  //cummulative input equation
};

/**
 * @brief Tree walker with information about edges source and target geometry
 * 
 * This class extens TargetWalker with information about the edges source node and offers an interface
 * to access the stored primitive geometry.
 * @tparam Kernel The math Kernel in use
 * @tparam SourcePrimitive The primitive geometry type stored at the edges source vertex
 * @tparam TargetPrimitive The primitive geometry type stored at the edges target vertex
 */
template<typename Kernel, template<class> class SourcePrimitive, template<class> class TargetPrimitive>
struct SourceTargetWalker : public TargetWalker<Kernel, TargetPrimitive>{
  
    SourceTargetWalker(const SourcePrimitive<Kernel>& sprim, const TargetPrimitive<Kernel>& tprim) 
                            : TargetWalker<Kernel, TargetPrimitive>(tprim), m_sourcePrimitive(sprim) {};
    
    /**
     * @brief Access the primitive geometry stored at the edges source vertex
     * @return const SourcePrimitive< Kernel >& The stored primitive geometry
     */
    const SourcePrimitive<Kernel>& getSourcePrimitive() {return m_sourcePrimitive;};

private:
    const SourcePrimitive<Kernel>&    m_sourcePrimitive; //primitive holding the value

};

/**
 * @brief Connects two nodes in a conditional manner
 *
 * This class describes the connection between two Nodes in the reduction tree. Connections are conditional,
 * they are only used as connections if some criteria is meet. To check this criteria the method 
 * \ref apply is used, it returns the validity of the connection. Furthermore an edge can have a own
 * behavior, an action can be executed. For this is custom class has to be derived which overrides 
 * the provided apply function.
 */
struct Connection {

    Node* source;
    Node* target;

    virtual ~Connection() {};
    
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

/**
 * @brief Exposes a Functor as connditional connection
 * 
 * This class provides a interface to use a functor as connection between nodes. The functor is called
 * and must behave like a normal connection: It takes the TreeWalker pointer as argument and returns
 * a boolean dependend if the criteria is met or not. The further traversal is handled by this class.
 * A possible implementation of a functor as lambda could look like this:
 * @code{.cpp}
 * auto functor = [](TreeWalker* walker) -> bool {return true;}; 
 * @endcode
 */
template<typename Functor>
struct FunctorConnection : public Connection {
  
    FunctorConnection(const Functor& f) : m_functor(f) {};
    
    virtual bool apply(TreeWalker* walker) const {
        
        if(m_functor(walker))
            return Connection::apply(walker);
        
        return false;
    };
    
private:
    const Functor& m_functor;
};

/**
 * @brief Helper class to easily create a constraint type and value checked connection
 * 
 * Often a connection between geometry nodes is conditioned on the type of the available constraint 
 * and the specific value of the constraint. To reduce the boilerplate needed for a connection of this 
 * common type this class provides a interface to define the allowed constrtaint type and the value it 
 * must have for the connection to be regarded as valid. This is specifed via the class template 
 * parameters. It can than be used with a functor which gets only supplied the primitive constraint to
 * do its action and which must not handle the validity check or the return type.
 * This class is intended to be used soley with Node::connectConditional. 
 * 
 * @note It is needed for the class to get the Final dcm system as template parameter. Hence it can 
 *       only be used in the context of a module 
 * 
 * A example of the inteded usage:
 * @code{.cpp}
 * template<typename Constraint, typename utilities::non_floating<typename Constraint::PimaryOptionType>::type option>
 * using EqualValue = ConstraintEqualValue<Final, Constraint, option>;
 * 
 * Node* mynode = ...
 * mynode->connectConditional<EqualValue<dcm::Distance, 0>>([](const dcm::Distance&) { ...;});
 * @endcode
 * 
 * @note The value cannot be a floating point number! If the option type is floating it gets convertet 
 *       to integer. This is a restriction coming from the c++ compiler which dies not allow floating 
 *       point numbers as type template arguments
 * 
 * @see ConstraintUnequalValue
 * @tparam Final The fully qualified dcm::system in use 
 * @tparam Constraint The primitive constraint type for which the connection should be used 
 * @tparam option The value of the primary optionn of the primitive constraint which must be met
 */
template<typename Final, typename Constraint, typename utilities::non_floating<typename Constraint::PimaryOptionType>::type option>
struct ConstraintEqualValue {
    
    template<typename Functor>
    struct Type : public Connection {
        
        Type(const Functor& f) : m_functor(f) {};
        
        virtual bool apply(TreeWalker* walker) const {
            
            auto cwalker = static_cast<ConstraintWalker<typename Final::Kernel>*>(walker);
            auto cons = cwalker->template getConstraint<Constraint>(Final::template constraintIndex<Constraint>::value);
            if(!cons)
                return false;
            
            if(cons->getPrimitiveConstraint().template getOption<0>() != option)
                return false;
            
            m_functor(cons->getPrimitiveConstraint());            
            cwalker->acceptConstraint(cons);
            return Connection::apply(walker);
        };
        
    private:
        const Functor& m_functor;
    };
};

/**
 * @brief Helper class to easily create a constraint type and value checked connection
 * 
 * This class works like ConstraintUnequalValue with the difference for searching for unequal constraint 
 * values. See \ref ConstraintEqualValue for documentation on usage.
 * 
 * @tparam Final The fully qualified dcm::system in use 
 * @tparam Constraint The primitive constraint type for which the connection should be used 
 * @tparam option The value of the primary optionn of the primitive constraint which must be met
 */
template<typename Final, typename Constraint, typename utilities::non_floating<typename Constraint::PimaryOptionType>::type option>
struct ConstraintUnequalValue {
    
    template<typename Functor>
    struct Type : public Connection {
        
        Type(const Functor& f) : m_functor(f) {};
        
        virtual bool apply(TreeWalker* walker) const {
            
            auto cwalker = static_cast<ConstraintWalker<typename Final::Kernel>*>(walker);
            auto cons = cwalker->template getConstraint<Constraint>(Final::template constraintIndex<Constraint>::value);
            if(!cons)
                return false;
            
            if(cons->template getOption<0>() == option)
                return false;
            
            m_functor(cons->getPrimitiveConstraint());                            
            cwalker->acceptConstraint(cons);
            return Connection::apply(walker);
        };
        
    private:
        Functor m_functor;
    };
};
    

/**
 * @brief Node of a reduction tree
 * 
 * Base class of reduction tree nodes. This class is used to build up connditional connection 
 * trees together with Conection class. It represents the reachable states of reduction. The class is 
 * responsible for holding all connections that go from it to other nodes. When a TreeWalker reaches
 * a Node it gets recorded inside the walker by setting its current node. Than the walker if applied
 * to all available connections to find a possible route to a next Node.
 * This class provides the interface for conecting Nodes together as well as for traversing the 
 * reduction tree by calling apply
 */
struct Node {
  
    virtual ~Node() {    
        //we own the edges,  hence we are responsible for deleting them
        for (Connection* e : m_edges)
            delete e;
    };
    
    /**
     * @brief Connect this node to annother via a specific Connection
     *
     * This function can be used to connect a node to annother one. As the connections are always 
     * conditional the Connection evaluating the condition needs to be supplied.
     * \Note The reduction tree is a static tree, hence there is no identifier for the made connection
     * and also no way to disconnect again
     *
     * \param node The node we build a connection to
     * \param edge The GeometryConnection which evaluates the condition for the transition
     */   
    template<typename Functor>
    void connect(Node* node, Functor func);    
    
    /**
    * @brief Add predicate determined connection
    * 
    * To reduce the boilerplate needed to write type and value checking functors for connections this 
    * method provides a interface for predicate based connections. It uses predicate classes which 
    * decide if a edge is valid or not and only pass the needed data to the functor provided to this 
    * method. A predicate must have the following structure, which allows the type initialisation by 
    * Functor type:
    * @code{.cpp}
    * struct Predicate {    *
    *    template<typename Functor>
    *    struct Type : public Connection {
    *    
    *    Type(const Functor& f) : m_functor(f) {};
    *    
    *    virtual bool apply(TreeWalker* walker) const {
    *        //check if this is a valid connection and return false if not...            
    *        m_functor(walker->MyFunkyDataExtraction);                            
    *        return true;
    *    };
    *    
    * private:
    *    Functor m_functor;
    * };
    * @endcode
    * 
    * This than allows to use the predicate with arbitrary functors, as long as those acceppt the argument
    * the predicate passes to it.
    */
    template<typename Predicate, typename Functor> 
    void connectConditional(Node* node, Functor func) {
        typedef typename Predicate::template Type<Functor> ConnectionType;
        connect(node, static_cast<Connection*>(new ConnectionType(func)));
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

        walker->setCurrentNode(this);
        for (Connection* e : m_edges) {
            if (e->apply(walker)) {
                return true;
                break;
            }
        };
        return false;
    }; 

private:
    std::vector<Connection*> m_edges;
};

template<typename Functor>
void Node::connect(Node* node, Functor func) {

    auto edge = new FunctorConnection<Functor>(func);
    connect(node, static_cast<Connection*>(edge));
};

template<>
inline void Node::connect<Connection*>(Node* node, Connection* edge) {

    edge->source = this;
    edge->target  = node;
    m_edges.push_back(edge);
};

//we can create the apply function of Connection only now as node must be fullydefined...
inline bool Connection::apply(TreeWalker* walker) const {
    
    target->apply(walker);
    return true;
};

/**
 * @brief Node for Geometries
 * 
 * This node is used for reduction of constraints, it provides extra functions to build equations
 * for the Geometry it represents. As a Node in a constraint reduction tree represents a geometry,
 * some times dependend, sometime undependend, it can be used to create the numeric equations for the 
 * symbolic reduction result.The provided interface is created to allow for late equation creation. The 
 * numeric Geometry must not be created at traversal time, but as the node reached in traversal is stored 
 * in the TreeWalker one can always access it later and use the functons provided here to create the 
 * appropriate numeric equations.
 */
template<typename Kernel>
struct GeometryNode : public Node {
    
    typedef shedule::FlowGraph::Node  FlowNode;
    typedef std::shared_ptr<numeric::Calculatable<Kernel>>  Equation;

    /**
     * @brief Creates numeric equation
     * 
     * @param walker The TreeWalker holding all relevant data needed for creating the numeric equation
     * @return Equation Numeric equation representing the nodes geometry type
     */
    virtual Equation buildGeometryEquation(TreeWalker* walker) const = 0;

    /**
     * @brief Creates numeric equation inside a \ref FlowGraph
     * 
     * This method create the numeric equation rigth in the provided flow graph and ensures the most
     * efficient execution of a recalculation by minimizing the virtual calls.
     * 
     * @param walker The TreeWalker holding all relevant data needed for creating the numeric equation
     * @return Equation Numeric equation representing the nodes geometry type
     */
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

/**
 * @brief ...
 * 
 */
struct EdgeReductionTree {

    typedef std::tuple<symbolic::Constraint*,symbolic::Geometry*,symbolic::Geometry*> SymbolicTuple;
    typedef std::vector<SymbolicTuple>                                                SymbolicVector;

    virtual ~EdgeReductionTree() {};
    
    /** @brief Analyses the global edges and finds the best reduction result
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
                                         const SymbolicVector& vec) = 0;
};

/**
 * @brief ...
 * 
 */
template<typename Kernel, template<class> class SourceGeometry, template<class> class TargetGeometry>
struct GeometryEdgeReductionTree : public EdgeReductionTree {

    typedef typename Kernel::Scalar Scalar;
    typedef EdgeReductionTree::SymbolicVector SymbolicVector;

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
                                         const SymbolicVector& symbolics) override {

        dcm_assert(source != target);
        //dcm_assert(dynamic_cast<TypeGeometry<Kernel, SourceGeometry>*>(source) != NULL);
        //dcm_assert(dynamic_cast<TypeGeometry<Kernel, SourceGeometry>*>(target) != NULL);

        //get the primitive geometries
        const SourceGeometry<Kernel>& pg1 = static_cast<symbolic::TypeGeometry<SourceGeometry<Kernel>>*>(source)->getPrimitve();
        const TargetGeometry<Kernel>& pg2 = static_cast<symbolic::TypeGeometry<TargetGeometry<Kernel>>*>(target)->getPrimitve();

        //create a new treewalker and set it up. We need to make sure that the correct walker for the 
        //given reduction job is used
        ConstraintWalker<Kernel>* walker = new SourceTargetWalker<Kernel, SourceGeometry, TargetGeometry>(pg1, pg2);
        walker->setConstraintPool(symbolics);        
        walker->setInitialNode(getSourceNode());        
        
        //start the calculation
        getSourceNode()->apply(walker);
        return walker;
    };
    
protected:
    reduction::UndependendGeometryNode<Kernel, TargetGeometry>* m_sourceNode;
    std::unordered_map<std::type_index, reduction::Node*>       m_nodesMap;
};

} //reduction

namespace numeric {

/**
    * @brief ...
    * 
    */
template<typename Kernel>
struct EdgeBuilder : public numeric::EquationBuilder<Kernel> {
    
    //the general reduction result: how many parameters can be eliminated
    int parameterReduction() {
        return m_paramReduction;
    }

private:
    int m_paramReduction; // the amount of degrees of freedom reduced by this reduction
};

/**
 * @brief ...
 * 
 */
template<typename Kernel>
struct ConstraintBuilder : public EdgeBuilder<Kernel> {
    
    typedef shedule::FlowGraph::Node                                                  Node;
    typedef std::shared_ptr<numeric::Calculatable<Kernel>>                            CalcPtr;
    typedef std::tuple<symbolic::Constraint*,symbolic::Geometry*,symbolic::Geometry*> SymbolicTuple;
    typedef std::vector<SymbolicTuple>                                                SymbolicVector;   
    typedef boost::multi_array<numeric::ConstraintEquationGenerator<Kernel>*,3>       GeneratorArray;
    
    ConstraintBuilder(const graph::LocalVertex& source, const graph::LocalVertex& target,
                      reduction::ConstraintWalker<Kernel>* sw, reduction::ConstraintWalker<Kernel>* tw,
                      const SymbolicVector& defaultConstraints, const GeneratorArray& cg) 
        : m_sourceVertex(source), m_targetVertex(target), m_sourceWalker(sw), m_targetWalker(tw), 
          m_constraints(defaultConstraints), m_generatorArry(cg) {};
    
    virtual ~ConstraintBuilder() {
        //we own the walkers
        delete m_sourceWalker;
        delete m_targetWalker;
    }
        
    //this function is used to create a default numeric geometry for the given symbolic one.
    virtual CalcPtr createGeometry(const graph::LocalVertex& vertex) override {
               
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
    createGeometryNode(const graph::LocalVertex& vertex, shedule::FlowGraph& flowgraph) override {
        
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
    virtual CalcPtr createReducedGeometry(const graph::LocalVertex& vertex) override {
        
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
    createReducedGeometryNode(const graph::LocalVertex& vertex, shedule::FlowGraph& flowgraph) override {
        
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
    createEquations(CalcPtr g1, CalcPtr g2) override {
        
        std::vector<CalcPtr> equations;
        for(auto tuple : m_constraints)
            equations.push_back(m_generatorArry[std::get<1>(tuple)->getType()]
                                               [std::get<2>(tuple)->getType()]
                                               [std::get<0>(tuple)->getType()]->buildEquation(g1, g2, std::get<0>(tuple)));
        
        return equations;
    };
                       
    //this function is used to create a node in the calculation flow graph with all default
    //constraint Equations
    virtual std::pair<std::vector<CalcPtr>, Node>
    createEquationsNode(CalcPtr g1, CalcPtr g2, shedule::FlowGraph& flow) override {
    
        if(m_constraints.size() == 1) {
            auto node = m_generatorArry[std::get<1>(m_constraints[0])->getType()]
                                       [std::get<2>(m_constraints[0])->getType()]
                                       [std::get<0>(m_constraints[0])->getType()]->buildEquationNode(g1, g2, 
                                                                                                std::get<0>(m_constraints[0]),
                                                                                                flow);
            std::vector<CalcPtr> vec(1);
            vec.push_back(node.first);
            return std::make_pair(vec, node.second);
        }
                
        //if we have multiple constraints we need to call the virtual functions anyway
        auto equations = createEquations(g1,g2);
        return std::make_pair(equations, flow.newActionNode([=](const shedule::FlowGraph::ContinueMessage& m){
            for(auto cons : equations)
                cons->execute();
        }));
    };

    //this function is used to create all default constraint equations
    virtual std::vector<CalcPtr>
    createReducedEquations(const graph::LocalVertex& target, CalcPtr g1, CalcPtr g2) override {
        
        std::vector<CalcPtr> equations;
        reduction::ConstraintWalker<Kernel>* walker = (m_targetVertex == target) ? m_targetWalker : m_sourceWalker;
        for(auto tuple : *walker)
            equations.push_back(m_generatorArry[std::get<1>(tuple)->getType()]
                                               [std::get<2>(tuple)->getType()]
                                               [std::get<0>(tuple)->getType()]->buildEquation(g1, g2, std::get<0>(tuple)));
        
        return equations;
    };
                       
    //this function is used to create a node in the calculation flow graph with all remaining
    //constraint Equations
    virtual std::pair<std::vector<CalcPtr>, Node>
    createReducedEquationsNode(const graph::LocalVertex& target, CalcPtr g1, CalcPtr g2, shedule::FlowGraph& flow) override {
        
        reduction::ConstraintWalker<Kernel>* walker = (m_targetVertex == target) ? m_targetWalker : m_sourceWalker;
        if(walker->size() == 1) {
            auto node = m_generatorArry[std::get<1>(walker->front())->getType()]
                                       [std::get<2>(walker->front())->getType()]
                                       [std::get<0>(walker->front())->getType()]->buildEquationNode(g1, g2, 
                                                                                               std::get<0>(walker->front()),
                                                                                               flow);
            std::vector<CalcPtr> vec(1);
            vec.push_back(node.first);
            return std::make_pair(vec, node.second);   
        }
                
        //if we have multiple constraints we need to call the virtual functions anyway
        auto equations = createReducedEquations(target,g1,g2);
        return std::make_pair(equations, flow.newActionNode([=](const shedule::FlowGraph::ContinueMessage& m){
            for(auto cons : equations)
                cons->execute();
        }));
    };
    
private:
    graph::LocalVertex                   m_sourceVertex, m_targetVertex;
    reduction::ConstraintWalker<Kernel>* m_sourceWalker;  //reduction result for the source geometry
    reduction::ConstraintWalker<Kernel>* m_targetWalker;  //reduction result for the target geometry 
    SymbolicVector                       m_constraints;   //all constraints in the edge together with its geometries
    const GeneratorArray&                m_generatorArry; //array with ConstraintGenerators for this special geometry combination
};
} //numeric

namespace symbolic {
    
/**
 * @brief Converts symbolic constraint systems into numeric equations
 * 
 * This class is used to analyse symbolic constraint systems and generate approprite numeric equations
 * for it. It searchews for the best possible numeric implementation by combining geometries and 
 * constraints appropriately. Note that at the time of reduction it is not yet known if the symbolic
 * system really should be converted to generalized coordinates, hence the numeric system is not 
 * directly created. Instead the NumericConverter adds an EquationBuilder to the graph which can be 
 * used to create the generalized or cartesian coordinate equations dependend on other analysis results.
 * 
 * @note This class is very heavy to initiate, it should not be created on the stack
 */
#include <boost/type_traits/is_same.hpp>
template<typename Kernel, typename GeometryList, typename ConstraintList, typename Graph>
struct NumericConverter {
              
    NumericConverter() {     
        
        int size = mpl::size<GeometryList>::type::value;
        m_treeArray.resize(boost::extents[size][size]);
        
        //build up the default reduction nodes
        //mpl trickery to get a sequence counting from 0 to the size of stroage entries
        typedef mpl::range_c<int,0,
                mpl::size<GeometryList>::value> StorageRange;
        //now iterate that sequence so we can access all storage elements with knowing the position
        //we are at
        utilities::RecursiveSequenceApplyer<GeometryList, ReductionTreeCreator> r(m_treeArray);
        mpl::for_each<StorageRange>(r);
        
        //do all the same for the constraints
        int constraints = mpl::size<ConstraintList>::type::value;
        m_generatorArray.resize(boost::extents[size][size][constraints]);        
        utilities::RecursiveSequenceApplyer<GeometryList, ConstraintGeneratorCreator<ConstraintList>::template type> g(m_generatorArray);
        mpl::for_each<StorageRange>(g);
    };
    
    ~NumericConverter() {
        //delete all reduction nodes
        int size = std::pow(mpl::size<GeometryList>::value, 2);
        //for(int i=0; i<size; ++i) 
        //    delete m_treeArray(i);
    };
    
    /**
     * @brief ...
     * 
     */
    void setupEquationBuilder(std::shared_ptr<Graph> g, graph::LocalEdge edge) {
        
        //get the geometry used in this edge
        symbolic::Geometry* source = g->template getProperty<symbolic::GeometryProperty>(g->source(edge));
        symbolic::Geometry* target = g->template getProperty<symbolic::GeometryProperty>(g->target(edge));
        
        dcm_assert(source);
        dcm_assert(target);
        
        //get the two reduction trees for this geometry combination
        reduction::EdgeReductionTree* stTree = m_treeArray[source->getType()][target->getType()];
        reduction::EdgeReductionTree* tsTree = m_treeArray[target->getType()][source->getType()];
     
        //get all constraints and cluster geometries, is needed
        typedef typename Graph::global_edge_iterator iterator;
        std::pair<iterator, iterator> it = g->getGlobalEdges(edge);
        std::vector<std::tuple<symbolic::Constraint*,symbolic::Geometry*,symbolic::Geometry*>> symbolics;
        for (; it.first != it.second; ++it.first) {
            auto c = g->template getProperty<symbolic::ConstraintProperty>(*it.first);
            //in case we have a cluster we need to access the global edge vertices directly
            auto sg = g->template getProperty<symbolic::GeometryProperty>((*it.first).source);
            auto tg = g->template getProperty<symbolic::GeometryProperty>((*it.first).target);
             
            symbolics.push_back(std::make_tuple(c, sg, tg));
        }
        
        //calculate both results
        auto* stWalker = static_cast<reduction::ConstraintWalker<Kernel>*>(stTree->apply(source, target, symbolics));
        auto* tsWalker = static_cast<reduction::ConstraintWalker<Kernel>*>(tsTree->apply(target, source, symbolics));
        
        //build the reduction and set store it in the graph. Make sure any pointer already stored is
        //deleted properly, especially the walkers
        auto reduction = g->template getProperty<numeric::EquationBuilderProperty<Kernel>>(edge);
        if(reduction)
            delete reduction;
        
        reduction = new numeric::ConstraintBuilder<Kernel>(g->source(edge), g->target(edge), 
                                                           tsWalker, stWalker, symbolics, m_generatorArray);
        g->template setProperty<numeric::EquationBuilderProperty<Kernel>>(edge, reduction);
    };
    
private:
    boost::multi_array<reduction::EdgeReductionTree*,2>                 m_treeArray;
    boost::multi_array<numeric::ConstraintEquationGenerator<Kernel>*,3> m_generatorArray;
       
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
        
            int idx1 = utilities::index<GeometryList, t1>::value;
            int idx2 = utilities::index<GeometryList, t2>::value;
            
            m_treeArray[idx1][idx2] = node1;
            m_treeArray[idx2][idx1] = node2;
        };
    };
    
 
    template<typename ConstraitSequence>
    struct ConstraintGeneratorCreator {
    
        template<typename GeometrySequence>
        struct type {
            boost::multi_array<numeric::ConstraintEquationGenerator<Kernel>*,3>& generator;
            
            type(boost::multi_array<numeric::ConstraintEquationGenerator<Kernel>*,3>& r) 
                : generator(r) {};
                
            template<typename N1, typename N2>
            void operator()() {
            
                typedef typename mpl::at<GeometrySequence, N1>::type G1;
                typedef typename mpl::at<GeometrySequence, N2>::type G2;
                
                typedef mpl::range_c<int,0,
                    mpl::size<ConstraitSequence>::value> StorageRange;
                //now iterate that sequence so we can access all storage elements with knowing the position
                //we are at
                InnerLoop<G1, G2, N1::value, N2::value> functor(generator);
                mpl::for_each<StorageRange>(functor);
            };
            
            template<typename G1, typename G2, int n1, int n2>
            struct InnerLoop {
                
                boost::multi_array<numeric::ConstraintEquationGenerator<Kernel>*,3>& generator;
            
                InnerLoop(boost::multi_array<numeric::ConstraintEquationGenerator<Kernel>*,3>& r) : generator(r) {};
                
                template<typename T>
                void operator()(const T& t) {
                
                    generator[n1][n2][T::value] = new numeric::TypedConstraintEquationGenerator<Kernel, 
                                                        typename mpl::at<ConstraintList, T>::type, G1, G2>();
                                                                            
                    generator[n2][n1][T::value] = generator[n1][n2][T::value];
                };
            };
        };
    };
};

}//symbolic
}//dcm

#endif //DCM_REDUCTION_H

