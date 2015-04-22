/*
    openDCM, dimensional constraint manager
    Copyright (C) 2012  Stefan Troeger <stefantroeger@gmx.net>

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

#ifndef FILTERGRAPH_HPP
#define FILTERGRAPH_HPP

#include "accessgraph.hpp"

#include <boost/graph/filtered_graph.hpp>

namespace dcm {
    
namespace details {

template<typename Graph>
struct group_filter {
    
    group_filter() {};
    group_filter(std::shared_ptr<Graph> g, int gr) : graph(g), group(gr) {};

    template<typename It>
    bool operator()(const It it) const {
        return (graph->template getProperty<graph::Group>(it) == group);
    };
    
private:
    int group;
    std::shared_ptr<Graph> graph;
};
    
template<typename Graph>
struct create_filtered_graph {
    
    template<typename T1, typename T2, typename T3, typename T4, typename T5>
    using type = boost::filtered_graph<boost::adjacency_list<T1,T2,T3,T4,T5>, 
                                    group_filter<Graph>, group_filter<Graph>>;
};

} //details
namespace graph {
    
    /** @addtogroup Core
 * @{
 * */

/**
* @ingroup ClusterGraph
* @brief A graph that represents a filtered subset of a clustergraph
* 
* This graph exposes only the edges and vertices of a given graph which belong to a 
* certain group. This allows to expose subsets of a graph. The FilterGraph implements
* the \ref AccessGraph api and hence allows full access to the structure and proeprties
* of the filtered graph. However, it does not provide any structure changing or 
* clustering functions. 
* 
* This class can and should be used to provide access to a subset of a given graph and 
* perform standart boost algorithms on it. It is possible to use multiple filtered 
* graphs with the same underlying graph in parallel if they do not filter for the same
* group
* 
* @see ClusterGraph
*/    
template<typename Graph>
class FilterGraph : public AccessGraph<typename Graph::edgeprop,
                                       typename Graph::globaledgeprop,
                                       typename Graph::vertexprop, 
                                       typename Graph::clusterprop, 
                                       details::create_filtered_graph<Graph>::template type> {
    
    typedef details::group_filter<Graph> Filter;
    
    typedef AccessGraph<typename Graph::edgeprop,
                        typename Graph::globaledgeprop,
                        typename Graph::vertexprop, 
                        typename Graph::clusterprop, 
                        details::create_filtered_graph<Graph>::template type> Base;
        
public:
    FilterGraph(std::shared_ptr<Graph> g, int group) : Base(m_graph), m_cluster(g), m_group(group),
                                      m_graph(g->getDirectAccess(), Filter(g, group), Filter(g, group)) {};
    
   /**
    * @brief A predicate object which decides whihc clusters belong to this filtered graph
    * 
    */
    struct ClusterFilter {
    
        std::shared_ptr<Graph> m_graph;
        int                    m_group;
        
        ClusterFilter(std::shared_ptr<Graph> g, int gr) : m_graph(g), m_group(gr) {};
        
        bool operator()(const typename  std::iterator_traits<typename Graph::cluster_iterator>::value_type& v) {
            return m_graph->template getProperty<Group>(v.first) == m_group;
        }
    };
    
    typedef boost::filter_iterator<ClusterFilter, typename Graph::cluster_iterator> cluster_iterator;
    typedef boost::filter_iterator<ClusterFilter, typename Graph::const_cluster_iterator> const_cluster_iterator;
                                              
    /**
     * @brief Iterators for all subclusters
     *
     * A pair with two \ref cluster_iterator is returned which point to the first cluster and
     * to one after the last. This allows full iteration over all subclusters. Note that the iterator does 
     * point to full clustergraphs, not filtered ones of any kind.
     *
     * @return :pair< cluster_iterator, cluster_iterator >
     **/
    std::pair<cluster_iterator, cluster_iterator> clusters();

    /**
     * @brief const equivalent to \ref clusters()
     *
     * A pair with two \ref cluster_iterator is returned which point to the first cluster and
     * to one after the last. This allows full iteration over all subclusters. Note that the iterator does 
     * point to full clustergraphs, not filtered ones of any kind.
     * 
     * @return :pair< const_cluster_iterator, const_cluster_iterator >
     **/
    std::pair<const_cluster_iterator, const_cluster_iterator> clusters() const;

    /**
     * @brief The amount of all subclusters
     *
     * Returns the number of all subclusters in this filtered graph
     * @return :size_t
     **/
    std::size_t numClusters() const;
    
private:
    typename Base::Graph   m_graph;
    std::shared_ptr<Graph> m_cluster;
    int                    m_group;
};

template<typename Graph>
std::pair<typename FilterGraph<Graph>::cluster_iterator, typename FilterGraph<Graph>::cluster_iterator> 
FilterGraph<Graph>::clusters() {
    
    ClusterFilter p = ClusterFilter(m_cluster, m_group);
    auto range = m_cluster->clusters();
    return std::make_pair(boost::make_filter_iterator(p, range.first, range.second),
                          boost::make_filter_iterator(p, range.second, range.second));
};

template<typename Graph>
std::pair<typename FilterGraph<Graph>::const_cluster_iterator, 
          typename FilterGraph<Graph>::const_cluster_iterator> FilterGraph<Graph>::clusters() const {
    
    ClusterFilter p = ClusterFilter(m_cluster, m_group);
    auto range = m_cluster->clusters();
    return std::make_pair(boost::make_filter_iterator(p, range.first, range.second),
                          boost::make_filter_iterator(p, range.second, range.second));
};

template<typename Graph>
std::size_t FilterGraph<Graph>::numClusters() const {
    
    int c=0;
    auto it = clusters();
    for(; it.first != it.second; ++it)
        ++c;
    
    return c;       
};

//convinience function for easy filter graph creation
template<typename Graph>
std::shared_ptr<FilterGraph<Graph>> make_filter_graph(std::shared_ptr<Graph> g, int group) {
    
    return std::make_shared<FilterGraph<Graph>>(g, group);
};
/** @} */
    
} //graph
} //dcm

#endif // FILTERGRAPH_HPP