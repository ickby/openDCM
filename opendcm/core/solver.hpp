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

#include "clustergraph.hpp"
#include "filtergraph.hpp"
#include "geometry.hpp"
#include "scheduler.hpp"
#include "reduction.hpp"

#include <boost/graph/connected_components.hpp>
#include <boost/graph/undirected_dfs.hpp>
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
    
    while(!done) {
        //start with edge analysing
        auto fedges = g->template filterRange<typename Graph::edge_changed>(g->edges());
        shedule::for_each(fedges.first, fedges.second, [&](graph::LocalEdge& e) {
            
            c.setupEquationBuilder(g, e);
            g->acknowledgeEdgeChanges(e);
        });
  
        //reduce edges into cluster if possible
        //TODO: create cluster from rigid edges
        
        //go on with cycle analysis and reduce into cluster if possible
        //TODO: find all cycles and analyse those for rigidy.
        
        //TODO: If we have reduced clusters we need to reiterate on the edges, so done ==false
        done = true;
    }
        
    //return c;
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

template<typename Graph, typename Kernel>
struct FlowGraphExtractor : public boost::default_dfs_visitor {

    typedef typename Graph::LocalVertex         LocalVertex;
    typedef typename Graph::LocalEdge           LocalEdge;
    typedef typename shedule::FlowGraph::Node   Node;

    shedule::FlowGraph flow;
//    std::map<LocalVertex, numeric::GeomertyEquation<Kernel>*> geometryMap;
    int parameter = 0;
    int equations = 0;
       
    void tree_edge(LocalEdge u, const Graph& graph) {
/*
         auto result = graph.template getProperty<symbolic::ResultProperty>(u);
         
         LocalVertex source = graph.source(u);
         auto res = geometryMap.find(source);
         if(res == geometryMap.end()) {
         
             res = geometryMap.find(graph.target(u));
             
             //if we still have nothing we need to create a start vertex geometry
             if(res == geometryMap.end()) {
                 
                 auto geom = graph.template getProperty<symbolic::GeometryProperty>(source);
                 //geometryMap[source] = result->(geom);
             }
         }
  */       
         
    };

    //back edges are special. if we are in clusters, all backedge constraints depend on all clusters
    //in the cluster cycle
    void back_edge(LocalEdge u, const Graph& graph) {

    
    };
};

/**
 * @brief Builds the numeric representation of the constraint system
 * 
 * This function builds up a executable which represents the geometric constraint system of a cluster 
 * in a numeric equaltion system. Executing the returned value recalculates all equations dependend on
 * the set parameters in the graph
 */
template<typename Kernel, typename Graph>
std::shared_ptr<shedule::Executable> buildGraphNumericSystem(std::shared_ptr<Graph> g) {
    
    //to ensure that we do not forget anything or that we have processed things twice we create a 
    //vector to store what waas processed so far.
    typedef std::shared_ptr<numeric::Calculatable<Kernel>> CalcPtr;
    std::map<graph::LocalVertex, CalcPtr> vProcessed;
    
    //we build up the numeric system for this graph. This also means finding parts that can be solved 
    //indivudual or simply after the main part (one-connected components)
    
    /*
    //to find one connected components we search all verices with one edge only. From there we can 
    //follow that edge to the adjacent vertices until we find a edge with more than two connections. That is
    //where the leaf ends.
    tbb::task_group leafs;
    tbb::concurrent_queue<tbb::flow::graph_node*> parallel_solvable;
    tbb::concurrent_queue<tbb::flow::graph_node*> end_solvable;
    
    int leafcount = 1e3;
    auto iter = g->vertices();
    for(;iter.first != iter.second; ++iter.first, ++leafcount) {
    
        //start the leaf processing if a vertex only has a single edge
        if(g->outDegree(*iter.first) == 1) {
            
            //doing this in parallel is not error free: it may happen that certain leafes are not found 
            //for example in a Y topology: when thread one removes one arm and thread two the other, it 
            //is possible that both detect the middle node as having 3 out edges and stop at the same 
            //position. Then the third arm stays behind as leaf.
            leafs.run([iter, g, leafcount, &s]() {                
                
                //group all edges and vertices that belong to the leaf
                graph::LocalVertex source = *iter.first;
                do {
                    g->template setProperty<graph::Group>(source, leafcount);
                    auto edges = g->outEdges(source);
                    graph::LocalEdge e;
                    if(g->template getProperty<graph::Group>(g->target(*edges.first)) != leafcount) 
                        e = *edges.first;
                    else 
                        e = *(edges.first++);
                    
                    g->template setProperty<graph::Group>(e, leafcount);
                    source = g->target(e);
                }
                while(g->outDegree(source) == 2);       
                       
                //See what kind of leafe we have here. Possibilities:
                //1. last leaf vertex is a cluster, than we can calculate the leaf in parallel and
                //   only the leaf-> main trunk connection needs to be calculated at the end
                //2. last leaf vertex is a geometry, than the whole leaf needs to be calculated 
                //   at the end
                if(g->template getProperty<graph::Type>(source) == graph::Cluster) {
                    
                }
            });            
        }
    }*/
   
    //shedule::FlowGraph* fg = new shedule::FlowGraph();
    auto fg = std::make_shared<numeric::CalculatableSequentialVector<Kernel>>();
    //iterate over all edges and create the equations for them
    auto edges = g->edges();
    for(; edges.first != edges.second; ++edges.first) {    
        typedef typename numeric::EquationBuilderProperty<Kernel> prop;
        numeric::EquationBuilder<Kernel>* builder = g->template getProperty<prop>(*edges.first); 
        auto source = g->source(*edges.first);
        CalcPtr g1;
        auto it = vProcessed.find(source);
        if( it == vProcessed.end() ) {
            g1 = builder->createGeometry(source);
            fg->addExecutable(g1);
            vProcessed[source] = g1;
        }
        else 
            g1 = it->second;
        
        auto target = g->target(*edges.first);
        CalcPtr g2;
        it = vProcessed.find(target);
        if( it == vProcessed.end() ) {
            g2 = builder->createGeometry(target);
            fg->addExecutable(g2);
            vProcessed[target] = g2;
        }
        else 
            g2 = it->second;
        
        auto vec = builder->createEquations(g1, g2);
        fg->append(vec);        
    }
   
    //build the solver object        
    return std::make_shared<numeric::Dogleg<Kernel>>(fg);
};
    

/**
 * @brief Creates a optimal numeric representation of the graph
 * 
 * This function returns a executable which can solve the dimensional constraint problem numerically.
 * It ensures that the system is optimally processed with symbolic methods and hence that the resulting
 * numeric equation system is minimal. Excecute will trigger the solving.
 */
template<typename Final, typename Graph, typename Converter>
std::shared_ptr<shedule::Executable> createSolvableSystem(std::shared_ptr<Graph> g, Converter& c) {

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
    auto iter = g->clusters();
    auto s1 = std::make_shared<shedule::ParallelConcurrentVector>();
    shedule::for_each(iter.first, iter.second, 
        [&](typename std::iterator_traits<typename Graph::cluster_iterator>::value_type& sub) {
            auto fg = createSolvableSystem<Final>(sub.second, c);
            s1->addExecutable(fg);
        }
    );        

    //now identify all individual components and create a executable for each
    std::shared_ptr<shedule::Executable> s2;
    if(components > 1) {       
        auto vec = std::make_shared<shedule::ParallelConcurrentVector>();
        shedule::for_each(0, components, 1, [&](int i) {
            auto filter = graph::make_filter_graph(g, i);
            vec->addExecutable(buildGraphNumericSystem<typename Final::Kernel>(filter));        
        });
        s2 = vec;
    }
    else {
        s2 = buildGraphNumericSystem<typename Final::Kernel>(g);   
    }
    
    
    //to calculate everthing we need to make sure that all subclusters are calculated before the 
    //main graph is handled.
    if(s1->empty()) {
        s1.reset();
        return s2;
    }    
    return shedule::make_executable([s1, s2](){
        s1->execute();
        s2->execute();
    });
}

} //solver

} //dcm
#endif //DCM_MODULE_CORE_H
