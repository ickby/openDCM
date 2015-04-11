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

template<typename EdgeFilter, typename VertexFilter>
struct create_filtered_graph {
    
    template<typename T1, typename T2, typename T3, typename T4, typename T5>
    using type = boost::filtered_graph<boost::adjacency_list<T1,T2,T3,T4,T5>, 
                                    EdgeFilter, VertexFilter>;
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
     */    
template<typename Graph, int Group>
class FilterGraph : public AccessGraph<typename Graph::edgeprop,
                                       typename Graph::globaledgeprop,
                                       typename Graph::vertexprop, 
                                       typename Graph::clusterprop, 
                                       details::create_filtered_graph<                                       
                                            typename Graph::template in_group<Group>,
                                            typename Graph::template in_group<Group>>::template type> {
    
    typedef typename Graph::template in_group<Group> Filter;
    typedef AccessGraph<typename Graph::edgeprop,
                        typename Graph::globaledgeprop,
                        typename Graph::vertexprop, 
                        typename Graph::clusterprop, 
                        details::create_filtered_graph<Filter, Filter>::template type> Base;
        
public:
    FilterGraph(std::shared_ptr<Graph> g) : Base(m_graph),
                                      m_graph(g->getDirectAccess(), Filter(g), Filter(g)) {};
    
private:
    typename Base::Graph m_graph;
};


/** @} */
    
} //graph
} //dcm

#endif // FILTERGRAPH_HPP