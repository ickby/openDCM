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

#ifndef DCM_MAPPING_H
#define DCM_MAPPING_H

#include <Eigen/Core>

#include "reduction.hpp"
                

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
    typedef std::map<symbolic::Geometry*, CalcPtr>         CalcMap;
    
    
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
    createGeometryNode(const graph::LocalVertex& vertex, std::shared_ptr<shedule::FlowGraph> flowgraph) = 0;
    
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
    virtual CalcPtr createReducedGeometry(const graph::LocalVertex& vertex, CalcPtr g) = 0;
     
    
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
    createReducedGeometryNode(const graph::LocalVertex& vertex, CalcPtr g, std::shared_ptr<shedule::FlowGraph> flowgraph) = 0;
    

    /**
     * @brief Creates numeric equations for all available constraints
     * 
     * This function creates numeric equations for all available binary constraints, no matter what the 
     * reduction did result in. The created equations of type BinaryEquation should have the geometries 
     * as input, which are at the EquationHandlers edge source/target. The geometries at source and target vertex
     * of the edge should be used as source and target equation for this function, mixing those will
     * lead to faulty equation claculation! 
     * 
     * @param map The map from symbolic geometries to their numeric equations
     * @return std::vector< CalcPtr > The vector with BinaryEquations of all constraints 
     */    
    virtual std::vector<CalcPtr>
    createBinaryEquations(std::map<symbolic::Geometry*, CalcPtr>& map) = 0;
         
    /**
     * @brief Creates numeric equations for the given vertex
     * 
     * This function creates numeric equations for all available unary constraints at the given vertex,
     * no matter what the reduction did result in. The created equations of type UnaryEquation have the 
     * geometry as input, which is provided with this function. It is the callers responsibility to ensure
     * that the given geometry equation is the one that belongs to the given vertex.
     * 
     * @param map The map from symbolic geometries to their numeric equations
     * @return std::vector< CalcPtr > The vector with BinaryEquation of all constraints 
     */    
    virtual std::vector<CalcPtr>
    createUnaryEquations(const graph::LocalVertex& vertec, std::map<symbolic::Geometry*, CalcPtr>& map) = 0;

    /**
     * @brief Creates numeric equations for all available constraints directly in a \ref FlowGraph
     * 
     * This method does exactly the same as \ref createEquations but creates those equations in the 
     * provided \ref FlowGraph. It is ensured that the calculation of the created equtions in the 
     * graph node is as efficient as possible. The returned node is unconnected. To allow initialisation 
     * of the equations created they are provided too.
     * 
     * @param map The map from symbolic geometries to their numeric equations
     * @param flow The \ref FlowGraph the node should be created in
     * @return std::pair< std::vector< CalcPtr >, Node > The constraint BinaryEquation vector and the Node 
     *                                                   in the \ref FlowGraph
     */
    virtual std::pair<std::vector<CalcPtr>, Node>
    createBinaryEquationsNode(std::map<symbolic::Geometry*, CalcPtr>& map, std::shared_ptr<shedule::FlowGraph> flow) = 0;
    
        /**
     * @brief Creates numeric equations for the given vertex directly in a \ref FlowGraph
     * 
     * This function creates numeric equations for all available unary constraints at the given vertex,
     * no matter what the reduction did result in. The created equations of type UnaryEquation have the 
     * geometry as input, which is provided with this function. It is the callers responsibility to ensure
     * that the given geometry equation is the one that belongs to the given vertex.
     * 
     * @param map The map from symbolic geometries to their numeric equations
     * @return std::pair< std::vector< CalcPtr >, Node > The constraint UnaryEquations vector and the Node 
     *                                                   in the \ref FlowGraph
     */    
    virtual std::pair<std::vector<CalcPtr>, Node>
    createUnaryEquationsNode(const graph::LocalVertex& vertex, std::map<symbolic::Geometry*, CalcPtr>& map,
                             std::shared_ptr<shedule::FlowGraph> flow) = 0;

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
     * @param map The map from symbolic geometries to their numeric equations
     * @return std::vector< CalcPtr > The vector with BinaryEquation of all remaining constraints 
     */  
    virtual std::vector<CalcPtr>
    createReducedEquations(const graph::LocalVertex& target, std::map<symbolic::Geometry*, CalcPtr>& map) = 0;
                       
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
     * @param map The map from symbolic geometries to their numeric equations
     * @param flow The \ref FlowGraph the node should be created in
     * @return std::pair< std::vector< CalcPtr >, Node > The constraint BinaryEquation vector and the Node 
     *                                                   in the \ref FlowGraph
     */
    virtual std::pair<std::vector<CalcPtr>, Node>
    createReducedEquationsNode(const graph::LocalVertex& target, std::map<symbolic::Geometry*, CalcPtr>& map,
                               std::shared_ptr<shedule::FlowGraph> flow) = 0;
        
         
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
     * @brief Returns the equations created for the processed edge
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

protected:
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
    typedef std::vector<BinarySymbolicTuple>                                          BinarySymbolicVector;   
    typedef std::vector<UnarySymbolicTuple>                                           UnarySymbolicVector; 
    typedef boost::multi_array<numeric::UnaryConstraintEquationGenerator<Kernel>*,2>  UnaryGeneratorArray;
    typedef boost::multi_array<numeric::BinaryConstraintEquationGenerator<Kernel>*,3> BinaryGeneratorArray;
    
    ConstraintEquationHandler(const graph::LocalVertex& source, const graph::LocalVertex& target,
                      reduction::ConstraintWalker<Kernel>* sw, reduction::ConstraintWalker<Kernel>* tw,
                      const BinarySymbolicVector& defaultConstraints, const UnarySymbolicVector& svConstraints, 
                      const UnarySymbolicVector& tvConstraints, const BinaryGeneratorArray& cg, 
                      const UnaryGeneratorArray& ucg, const std::vector<int>& idToIndex) 
        : m_sourceVertex(source), m_targetVertex(target), m_sourceWalker(sw), m_targetWalker(tw), 
          m_constraints(defaultConstraints), m_sourceConstraints(svConstraints), m_targetConstraints(tvConstraints),
          m_binaryGeneratorArray(cg), m_unaryGeneratorArray(ucg), m_idToIndex(idToIndex) {
              
        //we need to setup the parameter reduction
        auto init = std::static_pointer_cast<reduction::GeometryNode<Kernel>>(sw->getInitialNode())->staticParameterCount();
        auto end  = std::static_pointer_cast<reduction::GeometryNode<Kernel>>(sw->getFinalNode())->staticParameterCount();
        EdgeEquationHandler<Kernel>::m_paramReduction = init - end;
    };
    
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
    createGeometryNode(const graph::LocalVertex& vertex, std::shared_ptr<shedule::FlowGraph> flowgraph) override {
        
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
    virtual CalcPtr createReducedGeometry(const graph::LocalVertex& vertex, CalcPtr g) override {
        
        if(vertex == m_sourceVertex) {
            //build the correct input equation
            m_sourceWalker->setInputEquation(g);
            m_sourceWalker->transformCurrentToFinalNodeInput();
            
            auto node = std::static_pointer_cast<reduction::GeometryNode<Kernel>>(m_sourceWalker->getFinalNode());
            auto eqn = node->buildGeometryEquation(m_sourceWalker);
            Base::setVertexEquation(eqn, vertex);
            return eqn;
        }
        else if(vertex == m_targetVertex)  {
            //build the correct input equation
            //build the correct input equation
            m_targetWalker->setInputEquation(g);
            m_targetWalker->transformCurrentToFinalNodeInput();
            
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
    createReducedGeometryNode(const graph::LocalVertex& vertex, CalcPtr g, std::shared_ptr<shedule::FlowGraph> flowgraph) override {
        
        if(vertex == m_sourceVertex) {
            //build the correct input equation
            m_sourceWalker->setInputEquation(g);
            m_sourceWalker->transformCurrentToFinalNodeInput();
            
            auto node = std::static_pointer_cast<reduction::GeometryNode<Kernel>>(m_sourceWalker->getFinalNode());
            auto eqn = node->buildGeometryEquationNode(m_sourceWalker, flowgraph);
            Base::setVertexEquation(eqn.first, vertex);
            return eqn;
        }
        else if(vertex == m_targetVertex)  {
            //build the correct input equation
            m_targetWalker->setInputEquation(g);
            m_targetWalker->transformCurrentToFinalNodeInput();
            
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
    createBinaryEquations(std::map<symbolic::Geometry*, CalcPtr>& map) override {
        
        std::vector<CalcPtr> equations;
        for(auto tuple : m_constraints) {
            auto g1 = map[std::get<1>(tuple)];
            auto g2 = map[std::get<2>(tuple)];
            dcm_assert(g1 && g2);
            equations.push_back(m_binaryGeneratorArray[m_idToIndex[std::get<1>(tuple)->getType()]]
                                                      [m_idToIndex[std::get<2>(tuple)->getType()]]
                                                      [m_idToIndex[std::get<0>(tuple)->getType()]]->buildEquation(g1, g2, std::get<0>(tuple)));
        }
        
        Base::setEdgeEquations(equations);
        return equations;
    };
    
    virtual std::vector<CalcPtr>
    createUnaryEquations(const graph::LocalVertex& vertex, std::map<symbolic::Geometry*, CalcPtr>& map) override {
        
        std::vector<CalcPtr> equations;
        UnarySymbolicVector  constraints;
        if(vertex == m_sourceVertex)
            constraints = m_sourceConstraints;
        else 
            constraints = m_targetConstraints;
        
        for(auto tuple : constraints) {
            auto Geometry = map[std::get<1>(tuple)];
            dcm_assert(Geometry);
            
            auto* gen = m_unaryGeneratorArray[m_idToIndex[std::get<1>(tuple)->getType()]]
                                             [m_idToIndex[std::get<0>(tuple)->getType()]];
            
            if(!gen->applyToEquation(Geometry, std::get<0>(tuple)))
                equations.push_back(gen->buildEquation(Geometry, std::get<0>(tuple)));
        }
        
        //Base::setVertexEquations(equations);
        return equations;
    };
                       
    //this function is used to create a node in the calculation flow graph with all default
    //constraint Equations
    virtual std::pair<std::vector<CalcPtr>, Node>
    createBinaryEquationsNode(std::map<symbolic::Geometry*, CalcPtr>& map, std::shared_ptr<shedule::FlowGraph> flow) override {
                
        if(m_constraints.size() == 1) {
            auto g1 = map[std::get<1>(m_constraints[0])];
            auto g2 = map[std::get<2>(m_constraints[0])];
            dcm_assert(g1 && g2);
            
            auto node = m_binaryGeneratorArray[m_idToIndex[std::get<1>(m_constraints[0])->getType()]]
                                       [m_idToIndex[std::get<2>(m_constraints[0])->getType()]]
                                       [m_idToIndex[std::get<0>(m_constraints[0])->getType()]]->buildEquationNode(g1, g2, 
                                                                                                std::get<0>(m_constraints[0]),
                                                                                                flow);
            std::vector<CalcPtr> vec = {node.first};
            Base::setEdgeEquations(vec);
            return std::make_pair(vec, node.second);
        }
                
        //if we have multiple constraints we need to call the virtual functions anyway
        auto equations = createBinaryEquations(map);
        Base::setEdgeEquations(equations);
        return std::make_pair(equations, flow->newActionNode([=](const shedule::FlowGraph::ContinueMessage& m){
            for(auto cons : equations)
                cons->execute();
        }));
    };

    virtual std::pair<std::vector<CalcPtr>, Node>
    createUnaryEquationsNode(const graph::LocalVertex& vertex, std::map<symbolic::Geometry*, CalcPtr>& map,
                             std::shared_ptr<shedule::FlowGraph> flow) override {
        
        UnarySymbolicVector  constraints;
        if(vertex == m_sourceVertex)
            constraints = m_sourceConstraints;
        else 
            constraints = m_targetConstraints;
        
        if(constraints.size() == 1) {
            
            auto Geometry = map[std::get<1>(constraints[0])];
            dcm_assert(Geometry);
            
            auto gen = m_unaryGeneratorArray[m_idToIndex[std::get<1>(constraints[0])->getType()]]
                                            [m_idToIndex[std::get<0>(constraints[0])->getType()]];
                                            
            if(!gen->applyToEquation(Geometry, std::get<0>(constraints[0]))) {
                auto node = gen->buildEquationNode(Geometry, std::get<0>(constraints[0]), flow);
                std::vector<CalcPtr> vec = {node.first};
                Base::setEdgeEquations(vec);
                return std::make_pair(vec, node.second);
            }
            else return std::make_pair(std::vector<CalcPtr>(), flow->newActionNode([](const shedule::FlowGraph::ContinueMessage& m){}));
        }
                
        //if we have multiple constraints we need to call the virtual functions anyway
        auto equations = createUnaryEquations(vertex, map);
        //Base::setEdgeEquations(equations);
        return std::make_pair(equations, flow->newActionNode([=](const shedule::FlowGraph::ContinueMessage& m){
            for(auto cons : equations)
                cons->execute();
        }));
    };
    
    //this function is used to create all default constraint equations
    virtual std::vector<CalcPtr>
    createReducedEquations(const graph::LocalVertex& target, std::map<symbolic::Geometry*, CalcPtr>& map) override {
        
        std::vector<CalcPtr> equations;
        reduction::ConstraintWalker<Kernel>* walker = (m_targetVertex == target) ? m_targetWalker : m_sourceWalker;
        for(auto tuple : walker->binaryConstraintPool()) {
            auto g1 = map[std::get<1>(tuple)];
            auto g2 = map[std::get<2>(tuple)];
            dcm_assert(g1 && g2);
            
            equations.push_back(m_binaryGeneratorArray[m_idToIndex[std::get<1>(tuple)->getType()]]
                                               [m_idToIndex[std::get<2>(tuple)->getType()]]
                                               [m_idToIndex[std::get<0>(tuple)->getType()]]->buildEquation(g1, g2, std::get<0>(tuple)));
        }
        
        //TODO: Not sure if this mapping works
        for(auto tuple : walker->unaryConstraintPool()) {
            auto geometry = map[std::get<1>(tuple)];
            dcm_assert(geometry);
            
            auto gen = m_unaryGeneratorArray[m_idToIndex[std::get<1>(tuple)->getType()]]
                                            [m_idToIndex[std::get<0>(tuple)->getType()]];
            if(!gen->applyToEquation(geometry, std::get<0>(tuple)))
                equations.push_back(gen->buildEquation(geometry, std::get<0>(tuple)));
        }
        
        Base::setEdgeEquations(equations);
        return equations;
    };
                       
    //this function is used to create a node in the calculation flow graph with all remaining
    //constraint Equations
    virtual std::pair<std::vector<CalcPtr>, Node>
    createReducedEquationsNode(const graph::LocalVertex& target, std::map<symbolic::Geometry*, CalcPtr>& map,
                               std::shared_ptr<shedule::FlowGraph> flow) override {
        
        reduction::ConstraintWalker<Kernel>* walker = (m_targetVertex == target) ? m_targetWalker : m_sourceWalker;
        if((walker->binaryConstraintPool().size() == 1) && walker->unaryConstraintPool().empty()) {
            auto front = walker->binaryConstraintPool().front();
            
            auto g1 = map[std::get<1>(front)];
            auto g2 = map[std::get<2>(front)];
            dcm_assert(g1 && g2);
            
            auto node = m_binaryGeneratorArray[m_idToIndex[std::get<1>(front)->getType()]]
                                              [m_idToIndex[std::get<2>(front)->getType()]]
                                              [m_idToIndex[std::get<0>(front)->getType()]]->buildEquationNode(g1, g2, 
                                                                                          std::get<0>(front),
                                                                                          flow);
            std::vector<CalcPtr> vec = {node.first};
            Base::setEdgeEquations(vec);
            return std::make_pair(vec, node.second);   
        }
                
        //if we have multiple constraints we need to call the virtual functions anyway
        auto equations = createReducedEquations(target, map);
        Base::setEdgeEquations(equations);
        return std::make_pair(equations, flow->newActionNode([=](const shedule::FlowGraph::ContinueMessage& m){
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
    const UnaryGeneratorArray&           m_unaryGeneratorArray; //array of constraint generators for unary constraitns
    const std::vector<int>&              m_idToIndex;           //mapping from type ids to index within the arrays
};
    
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
        
        //setup the mapping from Id to index
        mpl::for_each<GeometryList>(IndexCreator(m_idToIndex));
        mpl::for_each<ConstraintList>(IndexCreator(m_idToIndex));
        mpl::for_each<UnaryConstraintList>(IndexCreator(m_idToIndex));
        
        //setup the geometric reduction trees
        int size = mpl::size<GeometryList>::type::value;
        m_treeArray.resize(boost::extents[size][size]);       
        utilities::RecursiveSequenceApplyer<GeometryList, ReductionTreeCreator> r(m_treeArray, m_idToIndex);
        mpl::for_each<GeometryList>(r);
        
        //setup the constraint generators
        int constraints = mpl::size<ConstraintList>::type::value;
        m_binaryGeneratorArray.resize(boost::extents[size][size][constraints]);        
        utilities::RecursiveSequenceApplyer<GeometryList, ConstraintGeneratorCreator<ConstraintList>> g(m_binaryGeneratorArray, m_idToIndex);
        mpl::for_each<GeometryList>(g);
        
        //and single geometry constraint generators
        constraints = mpl::size<UnaryConstraintList>::type::value;
        m_unaryGeneratorArray.resize(boost::extents[size][constraints]);        
        utilities::RecursiveSequenceApplyer<UnaryConstraintList, UnaryConstraintGeneratorCreator> gs(m_unaryGeneratorArray, m_idToIndex);
        mpl::for_each<GeometryList>(gs);
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
        reduction::EdgeReductionGraph* stTree = m_treeArray[m_idToIndex[source->getType()]][m_idToIndex[target->getType()]];
        reduction::EdgeReductionGraph* tsTree = m_treeArray[m_idToIndex[target->getType()]][m_idToIndex[source->getType()]];
     
        //get all constraints and cluster geometries, they are needed
        typedef typename Graph::global_edge_iterator iterator;
        std::pair<iterator, iterator> it = g->getGlobalEdges(edge);
        std::vector<std::tuple<symbolic::Constraint*,symbolic::Geometry*,symbolic::Geometry*>> stSymbolics, tsSymbolics;
        for (; it.first != it.second; ++it.first) {
            auto c = g->template getProperty<symbolic::ConstraintProperty>(*it.first);
            
            //in case we have a cluster we need to access the global edge vertices in the correct graph            
            auto res = g->getLocalVertexGraph((*it.first).source); 
            dcm_assert(fusion::at_c<2>(res));
            auto sg = fusion::at_c<1>(res)->template getProperty<symbolic::GeometryProperty>(fusion::at_c<0>(res));
            
            res = g->getLocalVertexGraph((*it.first).target);
            dcm_assert(fusion::at_c<2>(res));
            auto tg = fusion::at_c<1>(res)->template getProperty<symbolic::GeometryProperty>(fusion::at_c<0>(res));
             
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
                                                           targetSymbolics, m_binaryGeneratorArray, m_unaryGeneratorArray,
                                                           m_idToIndex);
        g->template setProperty<numeric::EquationHandlerProperty<Kernel>>(edge, reduction);
    };
    
    reduction::EdgeReductionGraph* getReductionGraph(int source, int target) {
        return m_treeArray[m_idToIndex[source]][m_idToIndex[target]];
    };
    
private:
    std::vector<int>                                                            m_idToIndex;
    boost::multi_array<reduction::EdgeReductionGraph*,2>                        m_treeArray;
    boost::multi_array<numeric::UnaryConstraintEquationGenerator<Kernel>*,2>    m_unaryGeneratorArray;
    boost::multi_array<numeric::BinaryConstraintEquationGenerator<Kernel>*,3>   m_binaryGeneratorArray;
       
    struct IndexCreator {
        std::vector<int>& idToIndex;
        int               counter = -1;
        
        IndexCreator(std::vector<int>& map) : idToIndex(map) {};
        
        template<typename T>
        void operator()(const T& t) {
            int id = T::id();
            if(idToIndex.size()<=(id+1))
                idToIndex.resize(id+1);
            
            idToIndex[id] = ++counter;
        };
    };
    
    struct ReductionTreeCreator {
    
        const std::vector<int>& m_idToIndex;
        boost::multi_array<reduction::EdgeReductionGraph*,2>& m_treeArray;
        
        ReductionTreeCreator(boost::multi_array<reduction::EdgeReductionGraph*,2>& r, const std::vector<int>& idToIndex) 
            : m_treeArray(r), m_idToIndex(idToIndex) {};
            
        template<typename T1, typename T2>
        void operator()() {
                    
            int idx1 = m_idToIndex[T1::id()];
            int idx2 = m_idToIndex[T2::id()];
            
            auto node1 = new reduction::GeometryEdgeReductionGraph<Kernel, 
                                geometry::extractor<T1>::template primitive,
                                geometry::extractor<T2>::template primitive >();
            
            m_treeArray[idx1][idx2] = node1;
            
            //if we have the same indexes we would override the first pointer and hence create a memory leak
            if(idx1 != idx2) {
                auto node2 = new reduction::GeometryEdgeReductionGraph<Kernel, 
                                    geometry::extractor<T2>::template primitive,
                                    geometry::extractor<T1>::template primitive >();
                  
                m_treeArray[idx2][idx1] = node2;
            }
        };
    };
    
 
    template<typename ConstraitSequence>
    struct ConstraintGeneratorCreator {
    
        const std::vector<int>& m_idToIndex;
        boost::multi_array<numeric::BinaryConstraintEquationGenerator<Kernel>*,3>& generator;
        
        ConstraintGeneratorCreator(boost::multi_array<numeric::BinaryConstraintEquationGenerator<Kernel>*,3>& r,
                                   const std::vector<int>& idToIndex) 
            : generator(r), m_idToIndex(idToIndex) {};
            
        template<typename G1, typename G2>
        void operator()() {
        
            InnerLoop<G1, G2> functor(generator, m_idToIndex);
            mpl::for_each<ConstraitSequence>(functor);
        };
        
        template<typename G1, typename G2>
        struct InnerLoop {
            
            const std::vector<int>& m_idToIndex;
            boost::multi_array<numeric::BinaryConstraintEquationGenerator<Kernel>*,3>& generator;
        
            InnerLoop(boost::multi_array<numeric::BinaryConstraintEquationGenerator<Kernel>*,3>& r, 
                      const std::vector<int>& idToIndex) : generator(r), m_idToIndex(idToIndex) {};
            
            template<typename Constraint>
            void operator()(const Constraint& t) {
                
                size_t idx1 = m_idToIndex[G1::id()];
                size_t idx2 = m_idToIndex[G2::id()];
                size_t idxC = m_idToIndex[Constraint::id()];
                generator[idx1][idx2][idxC] = new numeric::TypedBinaryConstraintEquationGenerator<Kernel, 
                                                    Constraint, G1, G2>();
                                                                        
                generator[idx2][idx1][idxC] = generator[idx1][idx2][idxC];
            };
        };
    };
    
    struct UnaryConstraintGeneratorCreator {
    
        const std::vector<int>& m_idToIndex;
        boost::multi_array<numeric::UnaryConstraintEquationGenerator<Kernel>*,2>& generator;
        
        UnaryConstraintGeneratorCreator(boost::multi_array<numeric::UnaryConstraintEquationGenerator<Kernel>*,2>& r,
                                        const std::vector<int>& idToIndex) 
                : generator(r), m_idToIndex(idToIndex) {};
                
        template<typename G, typename Constraint>
        void operator()() {
            
            size_t idxG = m_idToIndex[G::id()];
            size_t idxC = m_idToIndex[Constraint::id()];
            generator[idxG][idxC] = new numeric::TypedUnaryConstraintEquationGenerator<Kernel, 
                                                Constraint, G>();
        };
    };
};

}//numeric
}//dcm

#endif //DCM_MAPPING_H

