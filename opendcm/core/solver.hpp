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
#include "analyse.hpp"
#include "scheduler.hpp"

#include <boost/graph/connected_components.hpp>
#include <boost/graph/undirected_dfs.hpp>
#include <boost/fusion/include/vector.hpp>
#include <boost/fusion/include/at.hpp>
#include <boost/multi_array.hpp>

namespace fusion = boost::fusion;

namespace dcm {
       
namespace symbolic {
    
    /**
     * @brief Reduces the Graph to the smallest possible system
     * 
     * 
     * This function analyses the systems constraints and geometries and tries to find possible 
     * reductions, this means rigid subsections which are moved to their own clusters or simple
     * replacements of constraints with dependend geometry.
     * 
     * It furthermore finds all disconnected components in the graph and assigns each vertex and
     * edege to the components they belong to via their group property.
     */    
template<typename Final, typename Graph>
int reduceGraph(std::shared_ptr<Graph> g, 
                 const boost::multi_array<symbolic::EdgeReductionTree<Final>*,2>& reduction) {
    
    //start with edge analysing
    auto fedges = g->template filterRange<typename Graph::edge_changed>(g->edges());
    shedule::for_each(fedges.first, fedges.second, [&](graph::LocalEdge& e) {
        
        symbolic::Geometry* g1 = g->template getProperty<symbolic::GeometryProperty>(g->source(e));
        symbolic::Geometry* g2 = g->template getProperty<symbolic::GeometryProperty>(g->target(e));
        symbolic::EdgeReductionTree<Final>* tree = reduction[g1->type][g2->type];
        dcm_assert(tree);
        tree->apply(g, e); 
        g->acknowledgeEdgeChanges(e);
    });
    
    g->initIndexMaps();
    graph::property_map<graph::Group, Graph, graph::LocalVertex> gmap(g);
    graph::property_map<graph::Index, Graph, graph::LocalVertex> imap(g);
    graph::property_map<graph::Color, Graph, graph::LocalVertex> cmap(g);
    int c = boost::connected_components(g->getDirectAccess(), 
                                                 gmap, boost::vertex_index_map(imap).color_map(cmap));
    
    //connected components only assigns vertices, lets also assign edges to groups
    auto edges = g->edges();
    tbb::parallel_for_each(edges.first, edges.second, [&](graph::LocalEdge& e) {
        g->template setProperty<graph::Group>(e, g->template getProperty<graph::Group>(g->source(e)));
    });
    
    return c;
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

template<typename Graph>
shedule::FlowGraph buildGraphNumericSystem(std::shared_ptr<Graph> g) {
    
    /*
    //we build up the numeric system for this graph. This also means finding parts that can be solved 
    //indivudual or simply after the main part (one-connected components)
    
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
    
    shedule::FlowGraph fg;
    
    //iterate over all edges and create the equations for them
    auto edges = g->edges();
    for(; edges.first != edges.second; ++edges.first) {
    
        //auto result = g->template getProperty<symbolic::ResultProperty<Kernel>>(*edges.first);
        
    }
    
    return fg;
};
    
template<typename Final, typename Graph>
shedule::Executable* createSolvableSystem(std::shared_ptr<Graph> g, 
                    const boost::multi_array<symbolic::EdgeReductionTree<Final>*,2>& reduction) {
    
    
    
    //simplify the graph as much as possible. This has to be done before the subcluste processing as it is
    //possible that subclusters are groupt into yet annother subcluster
    int components = symbolic::reduceGraph(g, reduction);   
    
    //accesses all subclusters and handle them to make sure they are properly calculated before we try to 
    //reduce the toplevel cluster. A subcluster has to be reduced and solved imediatly, in contrary to the
    //toplevel system. This is a one time action, there is no need to store the subcluster calculation task.
    auto iter = g->clusters();
    tbb::parallel_for_each(iter.first, iter.second, 
        [&](typename std::iterator_traits<typename Graph::cluster_iterator>::value_type& sub) {
            auto fg = createSolvableSystem(sub.second, reduction);
            fg->execute();
            delete fg;
        }
    );        
        
    shedule::ParallelVector* s = new shedule::ParallelVector();
    //now identify all ndividual components and create a executable for each
    for(int i=0; i<components; ++i) {
        auto filter = graph::make_filter_graph(g, i);
        s->add(buildGraphNumericSystem(filter));        
    }
    
    //everything has been processed, lets return the solvable and let the caller decide what happens next
    return s;
}

} //solver

} //dcm
#endif //DCM_MODULE_CORE_H
