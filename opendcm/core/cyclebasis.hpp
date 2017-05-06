// Copyright Sebastien Kramm - 2015
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)



#ifndef DCM_CYCLEBASIS_CORE_HPP
#define DCM_CYCLEBASIS_CORE_HPP

#include <vector>

#include "accessgraph.hpp"

#include <boost/graph/undirected_dfs.hpp>

namespace dcm {
namespace graph {

template<typename Graph>
struct CycleBasis {
        
    CycleBasis(std::shared_ptr<Graph> g) : m_graph(g) {};
    
    void findCycles() {

        std::vector<LocalVertex> v_source_vertex;
        CycleDetector ld(v_source_vertex);

        // vertex and edge color maps
        m_graph->initIndexMaps();
        graph::property_map<graph::Color, Graph, graph::LocalEdge>   ecmap(m_graph);
        graph::property_map<graph::Color, Graph, graph::LocalVertex> vcmap(m_graph);
        
        boost::undirected_dfs( m_graph->getDirectAccess(), boost::visitor(ld).vertex_color_map(vcmap).edge_color_map(ecmap) );

        if ( !ld.cycleDetected() ) {
            m_cycles = std::vector<std::vector<LocalVertex>>(); // return empty vector, no cycles found
            return;
        }

        // else, get the cycles.
        std::vector<std::vector<LocalVertex>> v_cycles;

        for ( const auto& vi: v_source_vertex ) {
            std::vector<std::vector<LocalVertex>> v_paths;
            std::vector<LocalVertex> newv ( 1, vi ); // start by one of the filed source vertex
            v_paths.push_back ( newv );
            explore ( vi, v_paths, v_cycles );   // call of recursive function
        }

        // post process 1: remove the paths that are identical but reversed
        std::vector<std::vector<LocalVertex>> v_cycles2 = RemoveOppositePairs ( v_cycles );

        // post process 2: remove twin paths
        std::vector<std::vector<LocalVertex>> v_cycles3 = RemoveIdentical ( v_cycles2 );

        // post process 3: remove non-chordless cycles
        m_cycles = RemoveNonChordless ( v_cycles3 );
    }
    
    bool hasCycles() {
        return !m_cycles.empty();
    };
    
    
    
private:
    std::shared_ptr<Graph> m_graph;
    std::vector<std::vector<LocalVertex>> m_cycles;
    
    /**
    Recursive function, explores edges connected to \c v1 until we find a cycle
    \warning Have to be sure there \b is a cycle, else infinite recursion !
    v1: the starting vertex we want to explore
    v_cycles: this is where we store the paths that have cycles
    */
    bool explore(const LocalVertex& v1, std::vector<std::vector<LocalVertex>>& vv_paths,
                std::vector<std::vector<LocalVertex>>& v_cycles, int depth = 0) {
        
        ++depth;
        static int max_depth = std::max ( depth, max_depth );
        assert ( vv_paths.size() >0 );

        typename Graph::local_out_edge_iterator ei, ei_end;
        boost::tie ( ei, ei_end ) = m_graph->outEdges(v1);
        size_t edge_idx = 0;

        std::vector<LocalVertex> src_path = vv_paths[vv_paths.size()-1];

        bool found = false;
        for ( ; ei != ei_end; ++ei, ++edge_idx ) {

            bool b = false;
            LocalVertex v2a = m_graph->source(*ei);
            LocalVertex v2b = m_graph->target(*ei);

            if ( v2b == v1 && v2a == src_path[0] ) // we just found the edge that we started on, so no need to finish the current iteration, just move on.
                continue;

            std::vector<LocalVertex> newv ( src_path );
            bool AddNode = true;
            if ( newv.size() > 1 )
                if ( newv[ newv.size()-2 ] == v2b )
                    AddNode = false;

            if ( AddNode ) {
                if ( std::find ( newv.cbegin(), newv.cend(), v2b ) != newv.cend() ) {
                    newv.push_back ( v2b );
                    v_cycles.push_back ( newv );
                    return true;
                } else {
                    newv.push_back ( v2b );        // else add'em and continue
                    vv_paths.push_back ( newv );
                    b = explore ( v2b, vv_paths, v_cycles, depth );
                }
            }
            if ( b )
                found = true;
        }
        return found;
    }

    /**
    Remove twins : vector that are the same, but in reverse order
    */
    std::vector<std::vector<LocalVertex>>
    RemoveOppositePairs ( const std::vector<std::vector<LocalVertex>>& v_cycles ) {
        assert ( v_cycles.size() );
        std::vector<std::vector<LocalVertex>> out;         // output vector
        std::vector<bool> flags ( v_cycles.size(), true ); // some flags to keep track of which elements are reversed

        for ( size_t i=0; i<v_cycles.size()-1; ++i ) {
            if ( flags[i] ) {
                std::vector<LocalVertex> rev = v_cycles[i];              // step 1: build a reversed copy of the current vector
                std::reverse ( rev.begin(), rev.end() );
                for ( size_t j=i+1; j<v_cycles.size(); ++j ) {           // step 2: parse the rest of the list, and check
                    if ( flags[j] && rev == v_cycles[j] ) {              // if similar, then
                        out.push_back ( v_cycles[i] );                   //  1 - add current vector into output
                        flags[j] = false;                                //  2 -  invalidate the reversed one
                    }
                }
            }
        }
        return out;
    }

    //-------------------------------------------------------------------------------------------
    void PutSmallestElemFirst ( std::vector<LocalVertex>& vec ) {
        auto it = std::min_element ( vec.begin(), vec.end() );    // rotate so that smallest is first
        std::rotate ( vec.begin(), it, vec.end() );
    }

    /**
    Helper function for RemoveIdentical()

    Given an input vector "DABCD", it will return "ABCD" (removal of duplicate element, and first element is the smallest)
    */
    template<typename T>
    std::vector<T>
    GetSortedTrimmed ( const std::vector<T>& v_in ) {
        assert ( v_in.front() == v_in.back() ); // check this is a cycle
        assert ( v_in.size() > 2 );           // a (complete) cycle needs to be at least 3 vertices long

        std::vector<T> v_out ( v_in.size() - 1 );                     // Trim: remove
        std::copy ( v_in.cbegin(), v_in.cend()-1, v_out.begin() );    // last element

        PutSmallestElemFirst ( v_out );

        if ( v_out.back() < v_out[1] ) {                  // if we have 1-4-3-2, then
            std::reverse ( v_out.begin(), v_out.end() );  // we transform it into 2-3-4-1
            PutSmallestElemFirst ( v_out );               // and put smallest first: 1-2-3-4
        }
        return v_out;
    }

    /**
    Remove identical strings that are the same up to the starting point
    It also sorts the paths by rotating them so that the node of smallest index is first
    */
    std::vector<std::vector<LocalVertex>>
    RemoveIdentical ( const std::vector<std::vector<LocalVertex>>& v_cycles ) {
        assert ( v_cycles.size() );

        if ( v_cycles.size() == 1 ) {                        // if single path in input, then we justs add it, after trimming/sorting
            std::vector<std::vector<LocalVertex>> out ( 1, GetSortedTrimmed ( v_cycles[0] ) );
            return out;
        }

        std::vector<std::vector<LocalVertex>> out ( v_cycles.size() );
        for ( size_t i=0; i<v_cycles.size(); i++ )           // 1 - fill output vector with sorted/trimmed paths
            out[i] = GetSortedTrimmed ( v_cycles[i] );

        std::sort ( out.begin(), out.end() );                // 2 - sort
        out.erase (                                          // 3 - erase the ones that are
            std::unique ( out.begin(), out.end() ),          //  consecutive duplicates
            out.end()
        );

        return out;
    }

    //-------------------------------------------------------------------------------------------
    /// Returns true if vertices \c v1 and \c v2 are connected
    /**
    http://www.boost.org/doc/libs/1_59_0/libs/graph/doc/IncidenceGraph.html#sec:out-edges
    */
    bool AreConnected ( const LocalVertex& v1, const LocalVertex& v2 ) {
        auto pair_edge = m_graph->outEdges(v1); // get iterator range on edges
        for ( auto it = pair_edge.first; it != pair_edge.second; ++it )
            if ( v2 == m_graph->target(*it) )
                return true;
        return false;
    }

    //-------------------------------------------------------------------------------------------
    /// Return true if cycle is chordless
    /**
    See
    - https://en.wikipedia.org/wiki/Cycle_(graph_theory)#Chordless_cycles

    Quote:
    "A chordless cycle in a graph, also called a hole or an induced cycle, is a cycle such that
    no two vertices of the cycle are connected by an edge that does not itself belong to the cycle."
    */
    bool IsChordless ( const std::vector<LocalVertex>& path ) {
        if ( path.size() < 4 ) // else, no need to test
            return true;

        for ( size_t i=0; i<path.size()-3; ++i ) {
            for ( size_t j=i+2; j<path.size()-1; ++j )

                if ( AreConnected ( path[i], path[j] ) )
                    return false;
        }
        return true;
    }

    //-------------------------------------------------------------------------------------------
    /// Third step, remove non-chordless cycles
    std::vector<std::vector<LocalVertex>>
    RemoveNonChordless ( const std::vector<std::vector<LocalVertex>>& v_in ) {
        std::vector<std::vector<LocalVertex>> v_out;
        v_out.reserve ( v_in.size() ); // to avoid unnecessary memory reallocations and copies
        for ( const auto& cycle: v_in )
            if ( IsChordless ( cycle ) )
                v_out.push_back ( cycle );
        return v_out;
    }

    //-------------------------------------------------------------------------------------------
    /// Cycle detector for an undirected graph
    /**
    Passed by value as visitor to \c boost::undirected_dfs()

    See http://www.boost.org/doc/libs/1_58_0/libs/graph/doc/undirected_dfs.html
    */
    struct CycleDetector : public boost::dfs_visitor<> {

        CycleDetector(std::vector<LocalVertex>& vec) :  v_source_vertex(vec) {}
        
        bool cycleDetected() const {
            return !v_source_vertex.empty();
        }
        
        void back_edge ( LocalEdge e, const typename Graph::Graph& g ) {  // is invoked on the back edges in the graph.
            LocalVertex vs = boost::source ( e, g );
            LocalVertex vt = boost::target ( e, g );

            if (                                                                                               // add vertex to
                std::find ( v_source_vertex.cbegin(), v_source_vertex.cend(), vs ) == v_source_vertex.cend()   // the starting point list
                &&                                                                                             // only if both are
                std::find ( v_source_vertex.cbegin(), v_source_vertex.cend(), vt ) == v_source_vertex.cend()   // not already inside
            )
                v_source_vertex.push_back ( vs );
        }
    private:
        std::vector<LocalVertex>& v_source_vertex;
    };
    
};

} //graph
} //dcm

#endif // DCM_CYCLEBASIS_CORE_HPP
