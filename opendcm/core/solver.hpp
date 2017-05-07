/*
    openDCM, dimensional constraint manager
    Copyright (C) 2015  Stefan Troeger <stefantroeger@gmx.net>

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

#ifndef DCM_SOLVER_CORE_H
#define DCM_SOLVER_CORE_H

#include "accessgraph.hpp"
#include "clustergraph.hpp"
#include "filtergraph.hpp"
#include "geometry.hpp"
#include "scheduler.hpp"
#include "reduction.hpp"
#include "cyclebasis.hpp"

#include <boost/graph/connected_components.hpp>
#include <boost/graph/undirected_dfs.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/fusion/include/vector.hpp>
#include <boost/fusion/include/at.hpp>
#include <boost/multi_array.hpp>

namespace fusion = boost::fusion;

namespace dcm {
       
namespace symbolic {
    
    /**
     * @brief Find all seperated graph components
     * 
     * This function finds all seperated graph parts and numbers them accordingly. This allows to 
     * treat the disconnected components individual and hence in parallel.
     * Creating a FilterGraph on the found component is possible by using incremental inidices 
     * starting from 0 till the returned component count is reached.
     * @note This function clears all groups, as only this way ensures that the components numbers are
     *       unique 
     * @note The indentified and numbered groups are all still in the same data structure. Hence 
     * modifing the structure of the unconnected components cannot be done in parallel but must be 
     * done sequentially. Hence it is sensible to only use this function after all structure changing 
     * algorithms have finished.
     */
    
template<typename Final, typename Graph>
int splitGraph(std::shared_ptr<Graph> g)   {
    
    g->clearGroups();
    g->initIndexMaps();
    graph::property_map<graph::Group, Graph, graph::LocalVertex> gmap(g);
    graph::property_map<graph::Index, Graph, graph::LocalVertex> imap(g);
    graph::property_map<graph::Color, Graph, graph::LocalVertex> cmap(g);
    int c = boost::connected_components(g->getDirectAccess(), 
                                                 gmap, boost::vertex_index_map(imap).color_map(cmap));
    
    //connected components only assigns vertices, lets also assign edges to groups
    //As we support multiple groups, one needs to assign the edges only to the groups 
    //shared by both vertices.
    auto edges = g->edges();
    shedule::for_each(edges.first, edges.second, [&](graph::LocalEdge& e) {
        std::vector<int> edge;
        auto v1 = g->template getProperty<graph::Group>(g->source(e));
        auto v2 = g->template getProperty<graph::Group>(g->target(e));
        std::sort(v1.begin(), v1.end());
        std::sort(v2.begin(), v2.end());

        std::set_intersection(v1.begin(),v1.end(),v2.begin(),v2.end(),std::back_inserter(edge));
        g->template setProperty<graph::Group>(e, edge);
    });
    
    return c;
};

/**
* @brief Reduces the Graph to the smallest possible system 
* 
* This function analyses the systems constraints and geometries and tries to find possible 
* reductions, this means rigid subsections which are moved to their own clusters or simple
* replacements of constraints with dependend geometry.
* 
*/    
template<typename , typename Graph, typename Converter>
void reduceGraph(std::shared_ptr<Graph> g, Converter& c) {
    
    bool done = false;
    graph::CycleBasis<Graph> basis(g);
    
    while(!done) {
        //start with edge analysing
        auto fedges = g->template filterRange<typename Graph::edge_changed>(g->edges());
        shedule::for_each(fedges.first, fedges.second, [&](graph::LocalEdge& e) {
            
            c.setupEquationHandler(g, e);
            g->acknowledgeEdgeChanges(e);
        });
  
        //reduce edges into cluster if possible
        //TODO: create cluster from rigid edges
        
        //go on with cycle analysis and reduce into cluster if possible
        basis.findCycles();
        //TODO: analyse cycles for rigidy.
        done = true; //TODO: If we have reduced clusters we need to reiterate on the edges, so done ==false
    }    
};

} //symbolic

namespace numeric {
  
    
template<typename Graph>
shedule::FlowGraph buildRecalculationFlow(std::shared_ptr<Graph> g) {
    
//     tbb::flow::continue_node< tbb::flow::continue_msg > fg;
//     return fg;
    
};
    
} //numeric
    
namespace solver {

struct Builder {
        
    template<typename AccessGraph, typename Kernel>
    struct FlowGraphExtractor : public boost::default_dfs_visitor {
        
        typedef std::shared_ptr<numeric::Calculatable<Kernel>> CalcPtr;
        typedef typename numeric::EquationHandlerProperty<Kernel> prop;

        FlowGraphExtractor( std::shared_ptr<shedule::FlowGraph> flow, 
                            std::shared_ptr<numeric::CalculatableSequentialVector<Kernel>> eqns,
                            std::vector<graph::LocalEdge> tree,
                            std::map<graph::LocalVertex, std::pair<CalcPtr, shedule::FlowGraph::Node>>& pMap,
                            std::shared_ptr<AccessGraph> graph) 
                : m_flow(flow), m_eqns(eqns), m_tree(tree), m_graph(graph), m_processMap(pMap) {};

        void run(graph::LocalVertex start) {
            
            startVertex(start);
            handleVertex(start);
        };
        
    private:
        void startVertex(graph::LocalVertex v) {

            auto e = *m_graph->outEdges(v).first;
            numeric::EquationHandler<Kernel>* builder = m_graph->template getProperty<prop>(e); 
            auto res = builder->createGeometryNode(v, m_flow);
            m_processMap[v] = res;
            
            auto eqn = builder->createUnaryEquationsNode(v, res.first, m_flow);
            if(!eqn.first.empty())
                m_flow->connect(res.second, eqn.second);
            
            //connect the start node 
            m_flow->connect(m_flow->getBroadcastNode(), res.second);
            
            //setup the eqns vector. must happen after unary creation, as it can change the parameter count
            m_eqns->addExecutable(res.first);
            m_eqns->append(eqn.first);
        };
        
        //we look at all out edges and process the ones that are in the tree and not yet handled
        void handleVertex(graph::LocalVertex v) {
            
            auto out = m_graph->outEdges(v);
            for(; out.first != out.second; ++out.first) {
             
                if(std::find(m_tree.begin(), m_tree.end(), *out.first) != m_tree.end()) {
                    m_tree.erase(std::remove(m_tree.begin(), m_tree.end(), *out.first), m_tree.end());
                    handleEdge(*out.first, v);
                }
            }
        };
                
        void handleEdge(graph::LocalEdge e, graph::LocalVertex start) {

            auto target = (m_graph->source(e)==start) ? m_graph->target(e) : m_graph->source(e);
            numeric::EquationHandler<Kernel>* builder = m_graph->template getProperty<prop>(e); 
            
            auto it = m_processMap.find(start);
            assert(it != m_processMap.end()); 
            std::pair<CalcPtr, shedule::FlowGraph::Node> g1 = it->second;
            
            it = m_processMap.find(target);
            assert( it == m_processMap.end() );
            std::pair<CalcPtr, shedule::FlowGraph::Node> g2 = builder->createReducedGeometryNode(target, g1.first, m_flow);
            m_processMap[target] = g2;
                      
            m_flow->connect(g1.second, g2.second);
            
            //create unary and binary equations
            auto ueqn = builder->createUnaryEquationsNode(target, g2.first, m_flow);
            if(!ueqn.first.empty())
                m_flow->connect(g2.second, ueqn.second);
            
            auto beqn = builder->createReducedEquationsNode(target, g1.first, g2.first, m_flow);
            if(!beqn.first.empty())
                m_flow->connect(g2.second, beqn.second);
            
            //setup the eqns vector. must happen after unary creation, as it can change the parameter count
            m_eqns->addExecutable(g2.first);
            m_eqns->append(ueqn.first);
            m_eqns->append(beqn.first);
            
            //go on iterating
            handleVertex(target);
        };
        
    private:
        std::shared_ptr<AccessGraph> m_graph;
        std::shared_ptr<shedule::FlowGraph> m_flow;
        std::shared_ptr<numeric::CalculatableSequentialVector<Kernel>> m_eqns;
        std::vector<graph::LocalEdge> m_tree;
        std::map<graph::LocalVertex, std::pair<CalcPtr, shedule::FlowGraph::Node>>& m_processMap;
    };
    
    class adress_writer {
    public:
        template <class VertexOrEdge>
        void operator()(std::ostream& out, const VertexOrEdge& v) const {
            out << "[label=\"" << std::internal << std::hex << std::setfill('0') << v << "\"]";
        }
    };

    /**
    * @brief Builds the numeric representation of the constraint system
    * 
    * This function builds up a executable which represents the geometric constraint system of a cluster 
    * in a numeric equaltion system. All possible optimisation must have been done befor this function
    * is called as no such thing happens here. This function creates the most optimal numeric system
    * to solve the given graph in a single action.
    * @note All Equation handler need to be setup corectly, they are used within this function.
    */
    template<typename Kernel, typename Graph>
    static void solveGraphNumericSystem(std::shared_ptr<Graph> g) {
        
        //maybe we don't need to do anything for example when the graph is a disconnected geometry or cluster
        if(g->edgeCount() == 0)
            return;
        
        //create the vertex index
        g->initIndexMaps();
        graph::property_map<graph::Index, Graph, graph::LocalVertex> imap(g);
        /*
        std::cout<<std::endl<<std::endl;
        boost::write_graphviz(std::cout, g->getDirectAccess(), adress_writer(), 
                              boost::default_writer(), boost::default_writer(), imap);
        std::cout<<std::endl<<std::endl;
        */
        //we build up the numeric system for this graph. This means finding the reduction sequence that
        //gives the minimal possible parameter count. For this we use a minimum spanning tree.
        
        //we need a weight map. Our parameter reduction comes in handy
        typedef std::map<graph::LocalEdge, unsigned int> WeightMap;
        WeightMap weight;
        boost::associative_property_map<WeightMap> propmapWeight(weight);
        auto edges = g->edges();
        for(;edges.first != edges.second; ++edges.first) {
            auto handler = g->template getProperty<numeric::EquationHandlerProperty<Kernel>>(*edges.first);
            weight[*edges.first] = 100-static_cast<numeric::EdgeEquationHandler<Kernel>*>(handler)->parameterReduction();
        }
              
        //and build the tree
        std::vector<graph::LocalEdge> spanningTree;
        boost::kruskal_minimum_spanning_tree(g->getDirectAccess(), 
                                             std::back_inserter(spanningTree), 
                                             boost::vertex_index_map(imap).weight_map(propmapWeight));


        //build the flowgraph from the spanning tree. 
        auto flow = std::make_shared<shedule::FlowGraph>();
        auto eqns = std::make_shared<numeric::CalculatableSequentialVector<Kernel>>();
        typedef std::shared_ptr<numeric::Calculatable<Kernel>> CalcPtr;
        std::map<graph::LocalVertex, std::pair<CalcPtr, shedule::FlowGraph::Node>> processMap;

        FlowGraphExtractor<Graph, Kernel> extractor(flow, eqns, spanningTree, processMap, g);
        extractor.run(*g->vertices().first);  //TODO: don't use starting point random as now but use fixed
        
        //the non spanning tree edges have not been build. let's do that now
        edges = g->edges();
        for(;edges.first != edges.second; ++edges.first) {
            
            auto edge = *edges.first;
            if(std::find(spanningTree.begin(), spanningTree.end(), edge) == spanningTree.end()) {
                typedef typename numeric::EquationHandlerProperty<Kernel> prop;
                numeric::EquationHandler<Kernel>* builder = g->template getProperty<prop>(edge); 
                
                auto g1 = processMap[g->source(edge)];
                auto g2 = processMap[g->target(edge)];
                
                auto res = builder->createBinaryEquationsNode(g1.first, g2.first, flow);
                if(!res.first.empty()) {
                    flow->connect(g1.second, res.second);
                    flow->connect(g2.second, res.second);
                }
                
                eqns->append(res.first);
            }
        }
            
        //build the solver object and solve
        numeric::Dogleg<Kernel> solver(eqns, flow);
        solver.calculate();
        
        //we need to make sure that after solving all results are written back
        edges = g->edges();
        for(; edges.first != edges.second; ++edges.first) {
              auto* builder = g->template getProperty<typename numeric::EquationHandlerProperty<Kernel>>(*edges.first); 
              builder->writeToSymbolic();
        }
    };
        

    /**
    * @brief Creates a optimal numeric representation of the graph
    * 
    * This function returns a executable which can solve the dimensional constraint problem numerically.
    * It ensures that the system is optimally processed with symbolic methods and hence that the resulting
    * numeric equation system is minimal. Excecute will trigger the solving.
    */
    template<typename Final, typename Graph, typename Converter>
    static void solveSystem(std::shared_ptr<Graph> g, Converter& c) {
       
        //first do all the needed preprocessing
        auto vit = g->vertices();
        std::for_each(vit.first, vit.second, [g](const graph::LocalVertex& v) {
            auto prop = g->template getProperty<details::GraphObjectProperty>(v);
            if(prop)
                prop->preprocessVertex(g, v, g->template getProperty<graph::VertexProperty>(v));
        });
                
        auto eit = g->edges();
        std::for_each(eit.first, eit.second, [g](const graph::LocalEdge& e) {
            auto geit = g->getGlobalEdges(e);
            std::for_each(geit.first, geit.second, [g](const graph::GlobalEdge& ge) {
                auto prop = g->template getProperty<details::GraphObjectProperty>(ge);
                if(prop)
                    prop->preprocessEdge(g, ge);
            });
        });
        
        //maybe we don't need to do anything. Preprocessing need to be done even if vertices only,
        //but from here we can omit anything when no constraints are available. (for example fix cluster)
        if(g->edgeCount() == 0)
            return;
        
        //simplify the graph as much as possible. This has to be done before the subcluste processing as it is
        //possible that subclusters are groupt into yet annother subcluster
        symbolic::reduceGraph<Final>(g, c);   
        //all topology changing actions are done, now we can find groups to split. TODO: Split can happen 
        //before reduceGraph and reduce can be done on every splitted subgraph. But currently FilterGraph 
        //cannot be used as it cant change the graph structure
        int components = symbolic::splitGraph<Final>(g);
        
        //accesses all subclusters and handle them to make sure they are properly calculated before we try to 
        //reduce the toplevel cluster. A subcluster has to be reduced and solved before the cluster graph is 
        //solved, as the graph depends on the finished subcluster values. However, all subclusters can be solved
        //in parallel
        //TODO: exclude fix cluster
        auto iter = g->clusters();
        shedule::for_each(iter.first, iter.second, 
            [&](typename std::iterator_traits<typename Graph::cluster_iterator>::value_type& sub) {
                solveSystem<Final>(sub.second, c);
            }
        );

        //now identify all individual disconnected components and solve them
        if(components>1) {
            shedule::for_each(0, components, 1, [&](int i) {
                    auto filter = graph::make_filter_graph(g, i);
                    solveGraphNumericSystem<typename Final::Kernel>(filter);        
            });
        }
        else 
            solveGraphNumericSystem<typename Final::Kernel>(g);
        
        //we cant postprocess now, as we don't know if we are the toplevel cluster        
    }
    
    template<typename Graph>
    static void postprocessSystem(std::shared_ptr<Graph> g) {
    
        //as we need to access abolutely every node in the cluste rstructure we do this recursive
        //first process all subcluster
        auto iter = g->clusters();
        std::for_each(iter.first, iter.second, 
            [&](typename std::iterator_traits<typename Graph::cluster_iterator>::value_type& sub) {
                postprocessSystem(sub.second);
            }
        );
        
        auto vit = g->vertices();
        std::for_each(vit.first, vit.second, [g](const graph::LocalVertex& v) {
            auto prop = g->template getProperty<details::GraphObjectProperty>(v);
            if(prop)
                prop->postprocessVertex(g, v, g->template getProperty<graph::VertexProperty>(v));
        });
        auto eit = g->edges();
        std::for_each(eit.first, eit.second, [g](const graph::LocalEdge& e) {
            auto geit = g->getGlobalEdges(e);
            std::for_each(geit.first, geit.second, [g](const graph::GlobalEdge& ge) {
                auto prop = g->template getProperty<details::GraphObjectProperty>(ge);
                if(prop)
                    prop->postprocessEdge(g, ge);
            });
        });
    };
};

} //solver

} //dcm
#endif //DCM_MODULE_CORE_H
