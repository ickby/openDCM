/*
    openGCM, geometric constraint manager
    Copyright (C) 2012  Stefan Troeger <stefantroeger@gmx.net>

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/

#ifndef GCM_PROPERTY_H
#define GCM_PROPERTY_H

#include <boost/graph/graph_traits.hpp>
#include <boost/fusion/sequence.hpp>
#include <boost/fusion/container/vector.hpp>

#include <boost/mpl/find.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/property_map/property_map.hpp>

namespace mpl = boost::mpl;
namespace fusion = boost::fusion;

namespace gcm {

struct vertex_property {};
struct edge_property {};
struct cluster_property {};

namespace details {

template<typename Graph>
struct vertex_selector {
    typedef typename boost::graph_traits<Graph>::vertex_descriptor key_type;
    typedef typename Graph::vertex_properties sequence_type;
};

template<typename Graph>
struct edge_selector {
    typedef typename boost::graph_traits<Graph>::edge_descriptor key_type;
    typedef typename Graph::edge_properties sequence_type;
};

template< typename Kind, typename Graph>
struct property_selector : public mpl::if_<boost::is_same<Kind, vertex_property>,
        vertex_selector<Graph>, edge_selector<Graph> >::type {};

template<typename T>
struct property_type {
    typedef typename T::type type;
};
template<typename T>
struct property_kind {
    typedef typename T::kind type;
};
}

template<typename T>
struct is_edge_property : boost::is_same<typename T::kind,edge_property> {};

template<typename T>
struct is_vertex_property : boost::is_same<typename T::kind,vertex_property> {};

template<typename T>
struct is_cluster_property : boost::is_same<typename T::kind,cluster_property> {};

template <typename Property, typename Graph>
class property_map  {

public:
    typedef typename gcm::details::property_selector<typename Property::kind, Graph>::key_type key_type;
    typedef typename Property::type value_type;
    typedef typename Property::type&  reference;
    typedef boost::lvalue_property_map_tag category;

    typedef Property property;
    typedef typename gcm::details::property_selector<typename Property::kind, Graph>::sequence_type sequence;

    property_map(Graph& g)
        : m_graph(g) { }

    Graph& m_graph;
};

template<typename P, typename G>
typename property_map<P,G>::value_type	get(const property_map<P,G>& map,
        typename property_map<P,G>::key_type key)  {

    typedef property_map<P,G> map_t;
    typedef typename mpl::find<typename map_t::sequence, typename map_t::property>::type iterator;
    typedef typename mpl::distance<typename mpl::begin<typename map_t::sequence>::type, iterator>::type distance;
    return  fusion::at<distance>(fusion::at_c<0>(map.m_graph[key]));
};

template <typename P, typename G>
void  put(const property_map<P,G>& map,
          typename property_map<P,G>::key_type key,
          const typename property_map<P,G>::value_type& value)  {

    typedef property_map<P,G> map_t;
    typedef typename mpl::find<typename map_t::sequence, typename map_t::property>::type iterator;
    typedef typename mpl::distance<typename mpl::begin<typename map_t::sequence>::type, iterator>::type distance;
    fusion::at<distance>(fusion::at_c<0>(map.m_graph[key])) = value;
};


template <typename P, typename G>
typename property_map<P,G>::reference at(const property_map<P,G>& map,
        typename property_map<P,G>::key_type key) {
    typedef property_map<P,G> map_t;
    typedef typename mpl::find<typename map_t::sequence, typename map_t::property>::type iterator;
    typedef typename mpl::distance<typename mpl::begin<typename map_t::sequence>::type, iterator>::type distance;
    return fusion::at<distance>(fusion::at_c<0>(map.m_graph[key]));
}


//now create some standart properties
//***********************************

struct type_prop {
      //states the type of a cluster
      typedef cluster_property kind;
      typedef int type;
};

struct changed_prop {
    typedef cluster_property kind;
    typedef bool type;
};

}


#endif //GCM_PROPERTY_H
