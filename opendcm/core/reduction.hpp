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

#include "defines.hpp"
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
 * This class is the entry point to the numeric solver setup, it is the only place were a symbolic graph 
 * node and its numeric representation meet, it is the interface between numeric and symbolic worlds. 
 * It is centered around the graph edges. As the setup of dependend geometry depends on the base as well
 * as the dependend geometry, having the numeric initializer at the edge comes natural. Only there, 
 * with information about source and target vertices as well as constraints, the equation building has 
 * all needed informations. 
 * It has all the functions needed to get the numeric geometries and constraints for a \ref ClusterGraph 
 * edge and its source/target vertices. As the usage of vertices as un- or dependend geometry is not yet 
 * decided at creation time, this class offers the possibility to create both, dependend on the callers 
 * needs. The caller is responsible for calling the appropriate functions for reduced/unreduced geometry
 * and constraints and to not mix them. Futhermore the class provides functions to create equations
 * directly in a \ref FlowGraph which allows for more efficient calculation than adding the equation to
 * a \ref FlowGraph after the creation. When using the provided methods a minimal amount of virtual calls 
 * is assured.
 * @note The returned equations are not yet initialised and hence cannot be used directly. It is the
 * callers responsibility to call \ref init for all equations
 * 
 * As this class is the only point where symbolic and numeric representations meet, it is not only 
 * required to create a equation from the graph, but also to write back all numeric results into the
 * graph. This means the created equations are getting stored, so that a mapping from vertex and edge
 * to equation is possible. Functions for accessing those equations by graph entity are provided. There 
 * is also a more general function available to do the backwriting \ref writeToSymbolic
 * @note The equations are stored as a unique entry for each graph entity, this means if you call the 
 *       creation methods multiple times they will be overriden and the older equations cannot be 
 *       mapped to the graph anymore
 * 
 */
template<typename Kernel>
struct EquationHandler {
    
    virtual ~EquationHandler(){};
    
    typedef shedule::FlowGraph::Node  Node;
    typedef std::shared_ptr<numeric::Calculatable<Kernel>> CalcPtr;
    
    
    /**
     * @brief Create a undependend geometry equation 
     * 
     * This function creates the undependend geometry equation for a given vertex. This vertex must 
     * be the source or target of the edge where the EquationHandler is added to. The result does 
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
     * EquationHandler is added to. The result does ignore all reduction results done on the edge. 
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
     * be the source or target of the edge where the EquationHandler is added to. The created equation 
     * is a \ref InputEquation depending on the edges other vertex geometry, which is not provided, and
     * the constraints which have been used for the reduction. Hence if this EquationHandler is
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
     * as input, which are at the EquationHandlers edge source/target. The geometries at source and target vertex
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
     * as input, which are at the EquationHandlers edge source/target. The geometries at source and 
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
        
         
    /**
     * @brief Writes the numeric quations back to their symbolic counterpart
     * When creating the numeric equations, the symbolic geometry values are writen to the numeric equations.
     * Once the solving is done,  the process needs to be reversed, the numeric values need to be passed to
     * the symbolic representations. This is done with this equation. It used the mapping between graph 
     * entity and equation that was build up during equation creation.
     */
    virtual void writeToSymbolic() = 0;
    
    /**
     * @brief Returns the equation created for given vertex
     * @param v Vertex descriptor for which to inquery the equation
     * @return CalcPtr The created equation, nullptr if not created yet
     */
    CalcPtr getVertexEquation(const graph::LocalVertex& v) {return m_vertexEquation[v];};
    
    /**
     * @brief Returns the equations created for th processed edge
     * 
     * @return std::vector< CalcPtr > Vector of equations, empty if non created yet
     */
    std::vector<CalcPtr> getEdgeEquations() {return m_edgeEquations;};
    
protected:
    void setVertexEquation(CalcPtr ptr, const graph::LocalVertex& v) {m_vertexEquation[v] = ptr;};
    void setEdgeEquations(std::vector<CalcPtr> v) {m_edgeEquations = v;};
    
private:
    std::map<graph::LocalVertex, CalcPtr> m_vertexEquation;
    std::vector<CalcPtr>                  m_edgeEquations;
};

/**
 * @brief Store a EquationHandler at a graph edge
 *
 * EquationHandler class is used to convert a symbolic edge with its symbolic vertices to numeric 
 * equations (and back). This builder is the result of an symbolic analysis and must be stored at 
 * the edge which was analysed. This property does make this possible, it is added as a LocalEdge 
 * property to the ClusterGraph. It stores the EquationHandler as pointer to allow polymoprphism. 
 * The initial value is a nullptr, hence always check for validity! 
 * 
 * \tparam Kernel The Kernel used 
 **/
template<typename Kernel>
struct EquationHandlerProperty {
    typedef EquationHandler<Kernel>* type;
    struct default_value {
        EquationHandler<Kernel>* operator()() {
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
   
/** @addtogroup Reduction
 * @{
 * */
    
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
     * @return std::shared_ptr<reduction::Node> Pointer to initial node
     */
    std::shared_ptr<reduction::Node> getInitialNode()      {return m_initialNode;};
    
    /** @brief Set initial node before starting the traversal
     */    
    void  setInitialNode(std::shared_ptr<reduction::Node> n) {m_initialNode = n;};

    /**
     * @brief Retrieve the node at which the tree traversal ended     * 
     * @return std::shared_ptr<reduction::Node> Pointer to final node
     */
    std::shared_ptr<reduction::Node> getFinalNode()      {return m_finalNode;};
    
    /**
     * @brief Set the node at which the walker currently is during traversal
     */
    void  setCurrentNode(std::shared_ptr<reduction::Node> n) {m_finalNode = n;};
    
private:
    std::shared_ptr<reduction::Node> m_initialNode;
    std::shared_ptr<reduction::Node> m_finalNode;
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
    typedef std::tuple<symbolic::Constraint*,symbolic::Geometry*,symbolic::Geometry*> BinarySymbolicTuple;
    typedef std::tuple<symbolic::Constraint*,symbolic::Geometry*>                     UnarySymbolicTuple;
    typedef std::vector<BinarySymbolicTuple>                                          BinarySymbolicVector;    
    typedef std::vector<UnarySymbolicTuple>                                           UnarySymbolicVector;  
    typedef typename BinarySymbolicVector::const_iterator                             Iterator;  
    typedef typename UnarySymbolicVector::const_iterator                              VertexIterator; 
    
    /** @brief Fill constraint pool with available symbolic constraints
     */
    void setConstraintPools(const BinarySymbolicVector& vec, const UnarySymbolicVector& vertexVec) {
        m_binaryConstraintPool = vec; 
        m_unaryConstraintPool = vertexVec;
    };
    
    //interface for accessing constraints
    const BinarySymbolicVector& binaryConstraintPool() {return m_binaryConstraintPool;};
    const UnarySymbolicVector&  unaryConstraintPool()  {return m_unaryConstraintPool;};
    
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
    symbolic::TypeConstraint<T>* getConstraint(int type, int arity) {
        if(arity == 2) {
            for(const auto& tuple : m_binaryConstraintPool) {
                if(std::get<0>(tuple)->getType() == type)
                    return static_cast<symbolic::TypeConstraint<T>*>(std::get<0>(tuple));
            }
        }
        else if(arity == 1) {
            //maybe this constraint exist at the vertex
            for(const auto& tuple : m_unaryConstraintPool) {
                if(std::get<0>(tuple)->getType() == type)
                    return static_cast<symbolic::TypeConstraint<T>*>(std::get<0>(tuple));
            }
        }
        return nullptr;
    };
    
    /**
     * @brief Get constraint tuple (constraint and geometries) by constraint type
     * 
     * @return const BinarySymbolicTuple& Tuple of constraint, source and target geometry
     */
    BinarySymbolicTuple getBinaryConstraintTuple(int ConstraintType) {
    
        for(auto tuple : m_binaryConstraintPool) {
            if(std::get<0>(tuple)->getType() == ConstraintType)
                return tuple;
        }
        return BinarySymbolicTuple();
    };
    
    /**
     * @brief Get constraint  vertex tuple (constraint and geometry) by constraint type
     * 
     * @return const UnarySymbolicTuple& Tuple of constraint and target geometry
     */
    UnarySymbolicTuple getUnaryConstraintTuple(int ConstraintType) {
    
        for(auto tuple : m_unaryConstraintPool) {
            if(std::get<0>(tuple)->getType() == ConstraintType)
                return tuple;
        }
        return UnarySymbolicTuple();
    };
    
    /**
     * @brief Accept constraint during tree traversal
     * 
     * When a constraint is used by an edge it can be consumed and than not be available anymore in
     * the constraint pool. This method is used for this: when called the given constraint is removed
     * from the pool
     */   
    void acceptConstraint(symbolic::Constraint* c) {
        //remove the constraint
        if(c->getArity() == 2) {
            for(const auto& tuple : m_binaryConstraintPool) {
                if(std::get<0>(tuple) == c) {
                    acceptConstraintTuple(tuple);
                    return;
                }
            }
        }
        else {
            //maybe this constraint exist at the vertex
            for(const auto& tuple : m_unaryConstraintPool) {
                if(std::get<0>(tuple) == c) {
                    acceptConstraintTuple(tuple);
                    break;
                }
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
    void acceptConstraintTuple(const BinarySymbolicTuple& tuple) {
        m_binaryConstraintPool.erase(std::remove(m_binaryConstraintPool.begin(), m_binaryConstraintPool.end(), tuple), m_binaryConstraintPool.end());
    }
    
    void acceptConstraintTuple(const UnarySymbolicTuple& tuple) {
        m_unaryConstraintPool.erase(std::remove(m_unaryConstraintPool.begin(), m_unaryConstraintPool.end(), tuple), m_unaryConstraintPool.end());
    }
    
private:
    BinarySymbolicVector m_binaryConstraintPool;  //all the double geometry constraints we want to reduce
    UnarySymbolicVector  m_unaryConstraintPool;  //all the single geometry constraints we want to reduce
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
    
    TargetWalker(Primitive<Kernel>& prim) : m_targetPrimitive(prim) {};
 
    /**
     * @brief Access the primitive geometry stored at the edges target vertex
     * @return const Primitive< Kernel >& The stored primitive geometry
     */
    Primitive<Kernel>& getTargetPrimitive() {return m_targetPrimitive;};
        
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
    Primitive<Kernel>&    m_targetPrimitive; //primitive holding the value
    Equation              m_inputEqn;  //cummulative input equation
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
  
    SourceTargetWalker(SourcePrimitive<Kernel>& sprim, TargetPrimitive<Kernel>& tprim) 
                            : TargetWalker<Kernel, TargetPrimitive>(tprim), m_sourcePrimitive(sprim) {};
    
    /**
     * @brief Access the primitive geometry stored at the edges source vertex
     * @return const SourcePrimitive< Kernel >& The stored primitive geometry
     */
    SourcePrimitive<Kernel>& getSourcePrimitive() {return m_sourcePrimitive;};

private:
    SourcePrimitive<Kernel>&    m_sourcePrimitive; //primitive holding the value

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

    std::shared_ptr<reduction::Node> source;
    std::shared_ptr<reduction::Node> target;

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
 * std::shared_ptr<reduction::Node> mynode = ...
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
            auto cons = cwalker->template getConstraint<Constraint>(Final::template constraintIndex<Constraint>::value,
                                                                    Final::template constraintIndex<Constraint>::arity);
            if(!cons)
                return false;
            
            if(cons->getPrimitive().template getOption<0>() != option)
                return false;
            
            m_functor(cons->getPrimitive());            
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
            auto cons = cwalker->template getConstraint<Constraint>(Final::template constraintIndex<Constraint>::value,
                                                                    Final::template constraintIndex<Constraint>::arity);
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
struct Node : std::enable_shared_from_this<Node> {
  
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
    void connect(std::shared_ptr<reduction::Node> node, Functor func);    
    
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
    void connectConditional(std::shared_ptr<reduction::Node> node, Functor func) {
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

        walker->setCurrentNode(shared_from_this());
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
void Node::connect(std::shared_ptr<reduction::Node> node, Functor func) {

    auto edge = new FunctorConnection<Functor>(func);
    connect(node, static_cast<Connection*>(edge));
};

template<>
inline void Node::connect<Connection*>(std::shared_ptr<reduction::Node> node, Connection* edge) {

    edge->source = shared_from_this();
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
 * The node is also responsible for writing back the numeric results to the symbolic representation.
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
    
    /**
     * @brief Write the numeric equation to the symbolic graph representation
     * 
     * @param walker The TreeWalker holding all relevant data needed for accessing the symbolic world
     * @param eqn The numeric equation whichs result shall be transfered to the symbolic world
     */
    virtual void writeToSymbolic(TreeWalker* walker, Equation eqn) = 0;
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
    
    void writeToSymbolic(TreeWalker* walker, Equation eqn) override {
        TargetWalker<Kernel, G>* gwalker = static_cast<TargetWalker<Kernel, G>*>(walker);
        auto geom = std::static_pointer_cast<numeric::Geometry<Kernel, G>>(eqn);
        gwalker->getTargetPrimitive() = geom->output();
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
    
    virtual Equation buildGeometryEquation(TreeWalker* walker) const override {
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
    
    virtual std::pair< Equation, FlowNode > buildGeometryEquationNode(TreeWalker* walker, shedule::FlowGraph& flowgraph) const override {
        
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
    
    void writeToSymbolic(TreeWalker* walker, Equation eqn) override {
        auto* gwalker = static_cast<TargetWalker<Kernel, geometry::extractor<Output>::template primitive>*>(walker);
        auto geom = std::static_pointer_cast<DerivedG>(eqn);
        gwalker->getTargetPrimitive() = geom->output();
    }; 
};

/**
 * @brief Connected-Node graph to analyse geometry constraint combinations
 * This grpah is the data structure to handle Node and Connection datatypes in a graph like manner. 
 * It stores the entry point of the graph traversal, the source node. All other nodes can be stored 
 * in this class, this is however optional. The EdgeReductionGraph also allows to start the graph
 * traversal for given symbolic geometries and constraints and returns the used TreeWalker.
 */
struct EdgeReductionGraph {

    typedef std::tuple<symbolic::Constraint*,symbolic::Geometry*,symbolic::Geometry*> BinarySymbolicTuple;
    typedef std::vector<BinarySymbolicTuple>                                                BinarySymbolicVector;
    typedef std::tuple<symbolic::Constraint*, symbolic::Geometry*>                    UnarySymbolicTuple;
    typedef std::vector<UnarySymbolicTuple>                                          UnarySymbolicVector;

    virtual ~EdgeReductionGraph() {};
    
    /** @brief Analyses the given geometry/constraint pair and find the best reduction
     *
     * This function tries to find the best reduction for the given geometry/constraint combination.
     * The source and target parameters are the symbolic geometries to replace by reductions. The 
     * constraints are given view the vec parameter of BinarySymbolicVector type. Note that this vector does
     * not only contain the constraints, but for each constraint the symbolic geometires it connects. 
     * Normally the vec geometries are equal to source and target, however, that is not the case if 
     * one connects a cluster. Than source or target are a cluster geometry and the constraint connect
     * different primitives like points or lines.
     * \note The caller owns the returned TreeWalker pointer.
     * 
     * @remark The function is reentrant but is not safe to be called on the same data from multiple threads
     *
     * @param source The symbolic source geometry
     * @param target The symbolic target geometry
     * @param vec The vector of symbolic tuples describing the constraint type and the connected geometries
     */
    virtual reduction::TreeWalker* apply(symbolic::Geometry* source, symbolic::Geometry* target,
                                         const BinarySymbolicVector& vec, 
                                         const UnarySymbolicVector& targetSymbolics) = 0;
                                         

    /**
    * @brief Get the initial node of the tree
    * Returns the inital node where the whole reduction process starts.
    * @return std::shared_ptr< dcm::reduction::Node > The starting node for all reductions
    */
    std::shared_ptr<reduction::Node> sourceNode() {return m_sourceNode;};
    
    /**
     * @brief Replces the initial reduction node
     * 
     * As the source node gets auto created by the tree it may not be what you want. This function allows
     * to replace the current source node with a custom one. 
     * \note The connections of the replaced node will not be copied over
     */
    void replaceSourceNode(std::shared_ptr<Node> node) {m_sourceNode = node;};
    
    /**
     * @brief Gets or creates the tree node for the given dependend geometry type
     * 
     * It is important to access nodes after the creation time to add new connections to it. To ease
     * this process the EdgeReductionGraph offers a way to store created nodes an access them by type. 
     * With getTreeNode<MyType>() one acceses the previously created node, or if not created yet, 
     * creates a new one. 
     * \note It is not enforced to create nodes this way, one can easily mixed custom created nodes 
     *       with the ones created and stored by this function. 
     */    
    template<typename T>
    std::shared_ptr<reduction::Node> getTreeNode() { 
        auto iter = m_nodesMap.find(std::type_index(typeid(T)));
        if(iter != m_nodesMap.end())
            return iter->second;
        
        m_nodesMap[std::type_index(typeid(T))] = std::make_shared<DependendGeometryNode<T>>();
        return m_nodesMap[std::type_index(typeid(T))];
    };
    
protected:
    std::shared_ptr<reduction::Node> m_sourceNode;
    std::unordered_map<std::type_index, std::shared_ptr<reduction::Node>> m_nodesMap;
};

/**
 * @brief ...
 * 
 */
template<typename Kernel, template<class> class SourceGeometry, template<class> class TargetGeometry>
struct GeometryEdgeReductionGraph : public EdgeReductionGraph {

    typedef typename Kernel::Scalar Scalar;
    typedef EdgeReductionGraph::BinarySymbolicVector BinarySymbolicVector;
    typedef EdgeReductionGraph::UnarySymbolicVector UnarySymbolicVector;

    GeometryEdgeReductionGraph() {
        //create the tree source note, it must always be available
        m_sourceNode = std::make_shared<reduction::UndependendGeometryNode<Kernel, TargetGeometry>>();
    };
    
    virtual ~GeometryEdgeReductionGraph() {};
    
    virtual reduction::TreeWalker* apply(symbolic::Geometry* source, symbolic::Geometry* target,
                                         const BinarySymbolicVector& symbolics, 
                                         const UnarySymbolicVector& targetSymbolics) override {

        dcm_assert(source != target);
        //dcm_assert(dynamic_cast<TypeGeometry<Kernel, SourceGeometry>*>(source) != NULL);
        //dcm_assert(dynamic_cast<TypeGeometry<Kernel, SourceGeometry>*>(target) != NULL);

        //get the primitive geometries
        SourceGeometry<Kernel>& pg1 = static_cast<symbolic::TypeGeometry<SourceGeometry<Kernel>>*>(source)->getPrimitve();
        TargetGeometry<Kernel>& pg2 = static_cast<symbolic::TypeGeometry<TargetGeometry<Kernel>>*>(target)->getPrimitve();

        //create a new treewalker and set it up. We need to make sure that the correct walker for the 
        //given reduction job is used
        ConstraintWalker<Kernel>* walker = new SourceTargetWalker<Kernel, SourceGeometry, TargetGeometry>(pg1, pg2);
        walker->setConstraintPools(symbolics, targetSymbolics);        
        walker->setInitialNode(sourceNode());        
        
        //start the calculation
        sourceNode()->apply(walker);
        return walker;
    };
};

// @}

} //reduction

namespace numeric {

/**
    * @brief ...
    * 
    */
template<typename Kernel>
struct EdgeEquationHandler : public numeric::EquationHandler<Kernel> {
    
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
struct ConstraintEquationHandler : public EdgeEquationHandler<Kernel> {
    
    typedef EdgeEquationHandler<Kernel>                                               Base;
    typedef shedule::FlowGraph::Node                                                  Node;
    typedef std::shared_ptr<numeric::Calculatable<Kernel>>                            CalcPtr;
    typedef std::tuple<symbolic::Constraint*,symbolic::Geometry*,symbolic::Geometry*> BinarySymbolicTuple;
    typedef std::tuple<symbolic::Constraint*,symbolic::Geometry*>                     UnarySymbolicTuple;
    typedef std::vector<BinarySymbolicTuple>                                                BinarySymbolicVector;   
    typedef std::vector<UnarySymbolicTuple>                                          UnarySymbolicVector; 
    typedef boost::multi_array<numeric::UnaryConstraintEquationGenerator<Kernel>*,2>  UnaryGeneratorArray;
    typedef boost::multi_array<numeric::BinaryConstraintEquationGenerator<Kernel>*,3> BinaryGeneratorArray;
    
    ConstraintEquationHandler(const graph::LocalVertex& source, const graph::LocalVertex& target,
                      reduction::ConstraintWalker<Kernel>* sw, reduction::ConstraintWalker<Kernel>* tw,
                      const BinarySymbolicVector& defaultConstraints, const UnarySymbolicVector& svConstraints, 
                      const UnarySymbolicVector& tvConstraints, const BinaryGeneratorArray& cg) 
        : m_sourceVertex(source), m_targetVertex(target), m_sourceWalker(sw), m_targetWalker(tw), 
          m_constraints(defaultConstraints), m_sourceConstraints(svConstraints), m_targetConstraints(tvConstraints),
          m_binaryGeneratorArray(cg) {};
    
    virtual ~ConstraintEquationHandler() {
        //we own the walkers
        delete m_sourceWalker;
        delete m_targetWalker;
    }
        
    //this function is used to create a default numeric geometry for the given symbolic one.
    virtual CalcPtr createGeometry(const graph::LocalVertex& vertex) override {
               
        if(vertex == m_sourceVertex) {
            auto node = std::static_pointer_cast<reduction::GeometryNode<Kernel>>(m_sourceWalker->getInitialNode());
            auto eqn = node->buildGeometryEquation(m_sourceWalker);
            Base::setVertexEquation(eqn, vertex);
            return eqn;
        }
        else if(vertex == m_targetVertex)   {
            auto node = std::static_pointer_cast<reduction::GeometryNode<Kernel>>(m_targetWalker->getInitialNode());
            auto eqn = node->buildGeometryEquation(m_targetWalker);
            Base::setVertexEquation(eqn, vertex);
            return eqn;
        }
        else 
            throw solving_error() <<  boost::errinfo_errno(41) << error_message("No equations can be created for given vertex");
    };

    virtual std::pair< CalcPtr, Node>
    createGeometryNode(const graph::LocalVertex& vertex, shedule::FlowGraph& flowgraph) override {
        
        if(vertex == m_sourceVertex) {
            auto node = std::static_pointer_cast<reduction::GeometryNode<Kernel>>(m_sourceWalker->getInitialNode());
            auto eqn = node->buildGeometryEquationNode(m_sourceWalker, flowgraph);
            Base::setVertexEquation(eqn.first, vertex);
            return eqn;
        }
        else if(vertex == m_targetVertex) {
            auto node = std::static_pointer_cast<reduction::GeometryNode<Kernel>>(m_targetWalker->getInitialNode());
            auto eqn = node->buildGeometryEquationNode(m_targetWalker, flowgraph);
            Base::setVertexEquation(eqn.first, vertex);
            return eqn;
        }
        else 
            throw solving_error() <<  boost::errinfo_errno(41) << error_message("No equations can be created for given vertex");
    };
    
    //this function creates a reduced numeric geometry equation. For this the equation it depends on 
    //is required
    virtual CalcPtr createReducedGeometry(const graph::LocalVertex& vertex) override {
        
        if(vertex == m_sourceVertex) {
            auto node = std::static_pointer_cast<reduction::GeometryNode<Kernel>>(m_sourceWalker->getFinalNode());
            auto eqn = node->buildGeometryEquation(m_sourceWalker);
            Base::setVertexEquation(eqn, vertex);
            return eqn;
        }
        else if(vertex == m_targetVertex)  {
            auto node = std::static_pointer_cast<reduction::GeometryNode<Kernel>>(m_targetWalker->getFinalNode());
            auto eqn = node->buildGeometryEquation(m_targetWalker);
            Base::setVertexEquation(eqn, vertex);
            return eqn;
        }
        else 
            throw solving_error() <<  boost::errinfo_errno(41) << error_message("No equations can be created for given vertex");
    };
     
    //this function is used to create a reduced numeric geometry equation and a node for it in the 
    //corresponding flow graph. For this the equation it depends on is required
    virtual std::pair< CalcPtr, Node>
    createReducedGeometryNode(const graph::LocalVertex& vertex, shedule::FlowGraph& flowgraph) override {
        
        if(vertex == m_sourceVertex) {
            auto node = std::static_pointer_cast<reduction::GeometryNode<Kernel>>(m_sourceWalker->getFinalNode());
            auto eqn = node->buildGeometryEquationNode(m_sourceWalker, flowgraph);
            Base::setVertexEquation(eqn.first, vertex);
            return eqn;
        }
        else if(vertex == m_targetVertex)  {
            auto node = std::static_pointer_cast<reduction::GeometryNode<Kernel>>(m_targetWalker->getFinalNode());
            auto eqn = node->buildGeometryEquationNode(m_targetWalker, flowgraph);
            Base::setVertexEquation(eqn.first, vertex);
            return eqn;
        }
        else 
            throw solving_error() <<  boost::errinfo_errno(41) << error_message("No equations can be created for given vertex");
    };
    
    //this function is used to create all default constraint equations
    virtual std::vector<CalcPtr>
    createEquations(CalcPtr g1, CalcPtr g2) override {
        
        std::vector<CalcPtr> equations;
        for(auto tuple : m_constraints)
            equations.push_back(m_binaryGeneratorArray[std::get<1>(tuple)->getType()]
                                               [std::get<2>(tuple)->getType()]
                                               [std::get<0>(tuple)->getType()]->buildEquation(g1, g2, std::get<0>(tuple)));
        
        Base::setEdgeEquations(equations);
        return equations;
    };
                       
    //this function is used to create a node in the calculation flow graph with all default
    //constraint Equations
    virtual std::pair<std::vector<CalcPtr>, Node>
    createEquationsNode(CalcPtr g1, CalcPtr g2, shedule::FlowGraph& flow) override {
    
        if(m_constraints.size() == 1) {
            auto node = m_binaryGeneratorArray[std::get<1>(m_constraints[0])->getType()]
                                       [std::get<2>(m_constraints[0])->getType()]
                                       [std::get<0>(m_constraints[0])->getType()]->buildEquationNode(g1, g2, 
                                                                                                std::get<0>(m_constraints[0]),
                                                                                                flow);
            std::vector<CalcPtr> vec(1);
            vec.push_back(node.first);
            Base::setEdgeEquations(vec);
            return std::make_pair(vec, node.second);
        }
                
        //if we have multiple constraints we need to call the virtual functions anyway
        auto equations = createEquations(g1,g2);
        Base::setEdgeEquations(equations);
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
        for(auto tuple : walker->binaryConstraintPool())
            equations.push_back(m_binaryGeneratorArray[std::get<1>(tuple)->getType()]
                                               [std::get<2>(tuple)->getType()]
                                               [std::get<0>(tuple)->getType()]->buildEquation(g1, g2, std::get<0>(tuple)));
        
        Base::setEdgeEquations(equations);
        return equations;
    };
                       
    //this function is used to create a node in the calculation flow graph with all remaining
    //constraint Equations
    virtual std::pair<std::vector<CalcPtr>, Node>
    createReducedEquationsNode(const graph::LocalVertex& target, CalcPtr g1, CalcPtr g2, shedule::FlowGraph& flow) override {
        
        reduction::ConstraintWalker<Kernel>* walker = (m_targetVertex == target) ? m_targetWalker : m_sourceWalker;
        if(walker->binaryConstraintPool().size() == 1) {
            auto front = walker->binaryConstraintPool().front();
            auto node = m_binaryGeneratorArray[std::get<1>(front)->getType()]
                                       [std::get<2>(front)->getType()]
                                       [std::get<0>(front)->getType()]->buildEquationNode(g1, g2, 
                                                                                          std::get<0>(front),
                                                                                          flow);
            std::vector<CalcPtr> vec(1);
            vec.push_back(node.first);
            Base::setEdgeEquations(vec);
            return std::make_pair(vec, node.second);   
        }
                
        //if we have multiple constraints we need to call the virtual functions anyway
        auto equations = createReducedEquations(target,g1,g2);
        Base::setEdgeEquations(equations);
        return std::make_pair(equations, flow.newActionNode([=](const shedule::FlowGraph::ContinueMessage& m){
            for(auto cons : equations)
                cons->execute();
        }));
    };
    
    //this function must make sure the numeric results are writen to the symbolic graph entities
    void writeToSymbolic() override {
        //the constraints have no intteresting results, hence we only write back the vertex equations
        auto eqn = Base::getVertexEquation(m_sourceVertex);
        if(eqn) {
            auto node = std::static_pointer_cast<reduction::GeometryNode<Kernel>>(m_sourceWalker->getFinalNode());
            node->writeToSymbolic(m_sourceWalker, eqn);
        }
        eqn = Base::getVertexEquation(m_targetVertex);
        if(eqn) {
            auto node = std::static_pointer_cast<reduction::GeometryNode<Kernel>>(m_targetWalker->getFinalNode());
            node->writeToSymbolic(m_targetWalker, eqn);
        }
    };
    
private:
    graph::LocalVertex                   m_sourceVertex, m_targetVertex;
    reduction::ConstraintWalker<Kernel>* m_sourceWalker;        //reduction result for the source geometry
    reduction::ConstraintWalker<Kernel>* m_targetWalker;        //reduction result for the target geometry 
    BinarySymbolicVector                 m_constraints;         //all constraints in the edge
    UnarySymbolicVector                  m_sourceConstraints;   //all constraints at the source vertex
    UnarySymbolicVector                  m_targetConstraints;   //all constraints at the target vertex
    const BinaryGeneratorArray&          m_binaryGeneratorArray;//array with ConstraintGenerators for this special geometry combination
};
} //numeric

namespace symbolic {
    
//base class to allow retrieving the numeric converter without knowing its type at runtime
struct NumericConverterBase {};
    
/**
 * @brief Converts symbolic constraint systems into numeric equations
 * 
 * This class is used to analyse symbolic constraint systems and generate approprite numeric equations
 * for it. It searchews for the best possible numeric implementation by combining geometries and 
 * constraints appropriately. Note that at the time of reduction it is not yet known if the symbolic
 * system really should be converted to generalized coordinates, hence the numeric system is not 
 * directly created. Instead the NumericConverter adds an EquationHandler to the graph which can be 
 * used to create the generalized or cartesian coordinate equations dependend on other analysis results.
 * 
 * @note This class is very heavy to initiate, it should not be created on the stack
 */
#include <boost/type_traits/is_same.hpp>
template<typename Kernel, typename GeometryList, typename ConstraintList, 
         typename UnaryConstraintList, typename Graph>
struct NumericConverter : NumericConverterBase {
              
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
        m_binaryGeneratorArray.resize(boost::extents[size][size][constraints]);        
        utilities::RecursiveSequenceApplyer<GeometryList, ConstraintGeneratorCreator<ConstraintList>::template type> g(m_binaryGeneratorArray);
        mpl::for_each<StorageRange>(g);
        
        //and single geometry constraints
        constraints = mpl::size<UnaryConstraintList>::type::value;
        m_unaryGeneratorArray.resize(boost::extents[size][constraints]);        
        utilities::SequenceApplyer<GeometryList, UnaryConstraintGeneratorCreator<UnaryConstraintList>::template type> gs(m_unaryGeneratorArray);
        mpl::for_each<StorageRange>(gs);
    };
    
    ~NumericConverter() {
        //delete all reduction graphs
        for(int i=0; i<m_treeArray.num_elements(); ++i) 
            delete m_treeArray.data()[i];
        
        //delete all constraint creators. be carefull as we use generators twice!
        int geometries = mpl::size<GeometryList>::value;
        int constraints = mpl::size<ConstraintList>::value;
        for(int i=0; i<geometries; ++i) {
            for(int j=i; j<geometries; ++j) {
                for(int k=0; k<constraints; ++k) {
                    delete m_binaryGeneratorArray[i][j][k];
                }
            }
        }
        
        //delete all single geometry constraints
        for(int i=0; i<m_unaryGeneratorArray.num_elements(); ++i) 
            delete m_unaryGeneratorArray.data()[i];
        
    };
    
    /**
     * @brief ...
     * 
     */
    void setupEquationHandler(std::shared_ptr<Graph> g, graph::LocalEdge edge) {
        
        //get the geometry used in this edge
        symbolic::Geometry* source = g->template getProperty<symbolic::GeometryProperty>(g->source(edge));
        symbolic::Geometry* target = g->template getProperty<symbolic::GeometryProperty>(g->target(edge));
        
        dcm_assert(source);
        dcm_assert(target);
        
        //get the two reduction trees for this geometry combination
        reduction::EdgeReductionGraph* stTree = m_treeArray[source->getType()][target->getType()];
        reduction::EdgeReductionGraph* tsTree = m_treeArray[target->getType()][source->getType()];
     
        //get all constraints and cluster geometries, they are needed
        typedef typename Graph::global_edge_iterator iterator;
        std::pair<iterator, iterator> it = g->getGlobalEdges(edge);
        std::vector<std::tuple<symbolic::Constraint*,symbolic::Geometry*,symbolic::Geometry*>> stSymbolics, tsSymbolics;
        for (; it.first != it.second; ++it.first) {
            auto c = g->template getProperty<symbolic::ConstraintProperty>(*it.first);
            //in case we have a cluster we need to access the global edge vertices directly
            auto sg = g->template getProperty<symbolic::GeometryProperty>((*it.first).source);
            auto tg = g->template getProperty<symbolic::GeometryProperty>((*it.first).target);
             
            stSymbolics.push_back(std::make_tuple(c, sg, tg));
            tsSymbolics.push_back(std::make_tuple(c, tg, sg));
        }
        
        //get all vertex constraint symbolics
        std::vector<std::tuple<symbolic::Constraint*,symbolic::Geometry*>> sourceSymbolics, targetSymbolics;
        for(symbolic::Constraint* constraint : g->template getProperty<symbolic::ConstraintListProperty>(g->source(edge)))
            sourceSymbolics.push_back(std::make_tuple(constraint, source));
        for(symbolic::Constraint* constraint : g->template getProperty<symbolic::ConstraintListProperty>(g->target(edge)))
            targetSymbolics.push_back(std::make_tuple(constraint, target));
        
        //calculate both results
        auto* stWalker = static_cast<reduction::ConstraintWalker<Kernel>*>(stTree->apply(source, target, stSymbolics, targetSymbolics));
        auto* tsWalker = static_cast<reduction::ConstraintWalker<Kernel>*>(tsTree->apply(target, source, tsSymbolics, sourceSymbolics));
        
        //build the reduction and set store it in the graph. Make sure any pointer already stored is
        //deleted properly, especially the walkers
        auto reduction = g->template getProperty<numeric::EquationHandlerProperty<Kernel>>(edge);
        if(reduction)
            delete reduction;
        
        reduction = new numeric::ConstraintEquationHandler<Kernel>(g->source(edge), g->target(edge), 
                                                           tsWalker, stWalker, stSymbolics, sourceSymbolics,
                                                           targetSymbolics, m_binaryGeneratorArray);
        g->template setProperty<numeric::EquationHandlerProperty<Kernel>>(edge, reduction);
    };
    
    reduction::EdgeReductionGraph* getReductionGraph(int source, int target) {
        return m_treeArray[source][target];
    };
    
private:
    boost::multi_array<reduction::EdgeReductionGraph*,2>                        m_treeArray;
    boost::multi_array<numeric::UnaryConstraintEquationGenerator<Kernel>*,2>    m_unaryGeneratorArray;
    boost::multi_array<numeric::BinaryConstraintEquationGenerator<Kernel>*,3>   m_binaryGeneratorArray;
       
    template<typename Sequence>
    struct ReductionTreeCreator {
    
        boost::multi_array<reduction::EdgeReductionGraph*,2>& m_treeArray;
        
        ReductionTreeCreator(boost::multi_array<reduction::EdgeReductionGraph*,2>& r) 
            : m_treeArray(r) {};
            
        template<typename N1, typename N2>
        void operator()() {
        
            typedef typename mpl::at<Sequence, N1>::type t1;
            typedef typename mpl::at<Sequence, N2>::type t2;
            
            int idx1 = utilities::index<GeometryList, t1>::value;
            int idx2 = utilities::index<GeometryList, t2>::value;
            
            auto node1 = new reduction::GeometryEdgeReductionGraph<Kernel, 
                                geometry::extractor<t1>::template primitive,
                                geometry::extractor<t2>::template primitive >();
            
            m_treeArray[idx1][idx2] = node1;
            
            //if we have the same indexes we would override the first pointer and hence create a memory leak
            if(idx1 != idx2) {
                auto node2 = new reduction::GeometryEdgeReductionGraph<Kernel, 
                                    geometry::extractor<t2>::template primitive,
                                    geometry::extractor<t1>::template primitive >();
                  
                m_treeArray[idx2][idx1] = node2;
            }
        };
    };
    
 
    template<typename ConstraitSequence>
    struct ConstraintGeneratorCreator {
    
        template<typename GeometrySequence>
        struct type {
            boost::multi_array<numeric::BinaryConstraintEquationGenerator<Kernel>*,3>& generator;
            
            type(boost::multi_array<numeric::BinaryConstraintEquationGenerator<Kernel>*,3>& r) 
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
                
                boost::multi_array<numeric::BinaryConstraintEquationGenerator<Kernel>*,3>& generator;
            
                InnerLoop(boost::multi_array<numeric::BinaryConstraintEquationGenerator<Kernel>*,3>& r) : generator(r) {};
                
                template<typename T>
                void operator()(const T& t) {
                
                    generator[n1][n2][T::value] = new numeric::TypedBinaryConstraintEquationGenerator<Kernel, 
                                                        typename mpl::at<ConstraintList, T>::type, G1, G2>();
                                                                            
                    generator[n2][n1][T::value] = generator[n1][n2][T::value];
                };
            };
        };
    };
    
    template<typename ConstraitSequence>
    struct UnaryConstraintGeneratorCreator {
    
        template<typename GeometrySequence>
        struct type {
            boost::multi_array<numeric::UnaryConstraintEquationGenerator<Kernel>*,2>& generator;
            
            type(boost::multi_array<numeric::UnaryConstraintEquationGenerator<Kernel>*,2>& r) 
                : generator(r) {};
                
            template<typename N, int n>
            void operator()() {
            
                typedef typename mpl::at<GeometrySequence, N>::type G;
                
                typedef mpl::range_c<int,0,
                    mpl::size<ConstraitSequence>::value> StorageRange;
                //now iterate that sequence so we can access all storage elements with knowing the position
                //we are at
                InnerLoop<G, n> functor(generator);
                mpl::for_each<StorageRange>(functor);
            };
            
            template<typename G, int n>
            struct InnerLoop {
                
                boost::multi_array<numeric::UnaryConstraintEquationGenerator<Kernel>*,2>& generator;
            
                InnerLoop(boost::multi_array<numeric::UnaryConstraintEquationGenerator<Kernel>*,2>& r) : generator(r) {};
                
                template<typename T>
                void operator()(const T& t) {
                
                    generator[n][T::value] = new numeric::TypedUnaryConstraintEquationGenerator<Kernel, 
                                                        typename mpl::at<UnaryConstraintList, T>::type, G>();
                };
            };
        };
    };
};

}//symbolic
}//dcm

#endif //DCM_REDUCTION_H

