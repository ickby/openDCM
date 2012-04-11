/*
    <one line to give the program's name and a brief idea of what it does.>
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

#ifndef NS2_PROPERTY_H
#define NS2_PROPERTY_H

#include <boost/graph/graph_traits.hpp>
#include <boost/fusion/sequence.hpp>
#include <boost/fusion/container/vector.hpp>

#include <boost/mpl/find.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/property_map/property_map.hpp>

namespace mpl = boost::mpl;
namespace fusion = boost::fusion;

namespace dcm {
  
  struct vertex_property_tag {};
  struct edge_property_tag {};
  
namespace details {
  
  template< typename Kind, typename Graph>
  struct property_selector {
    typedef void key_type;
    typedef void bundle_type;
  };
  
  template<typename Graph>
  struct property_selector<vertex_property_tag, Graph> {
    typedef typename boost::graph_traits<Graph>::vertex_descriptor key_type;
    typedef typename boost::vertex_bundle_type<Graph>::type 	  bundle_type;
  };
  
  template<typename Graph>
  struct property_selector<edge_property_tag, Graph> {
    typedef typename boost::graph_traits<Graph>::edge_descriptor key_type;
    typedef typename boost::edge_bundle_type<Graph>::type 	bundle_type;
  };
}

template<typename T>
struct is_edge_property : boost::is_same<typename T::kind,edge_property_tag>{};

template<typename T>
struct is_vertex_property : boost::is_same<typename T::kind,vertex_property_tag>{};

  
template <typename Property, typename Graph>
class fusion_property_map  {
  
  public:
    typedef typename dcm::details::property_selector<typename Property::kind, Graph>::key_type key_type; 
    typedef typename Property::type value_type;
    typedef typename Property::type&  reference;
    typedef boost::lvalue_property_map_tag category;
    
    typedef Property property;
    typedef typename dcm::details::property_selector<typename Property::kind, Graph>::bundle_type bundle;

    fusion_property_map(Graph& g) 
      : m_graph(g) { }
    
    Graph& m_graph;
};
  
template<typename P, typename G>
typename fusion_property_map<P,G>::value_type	get(const fusion_property_map<P,G>& map,
						    typename fusion_property_map<P,G>::key_type key)  {
    
    typedef fusion_property_map<P,G> map_t;
    typedef typename mpl::find<typename map_t::bundle, typename map_t::property>::type iterator;
    return  fusion::at_c<iterator::pos::value>(fusion::at_c<0>(map.m_graph[key]));
};
  
template <typename P, typename G>
void  put(const fusion_property_map<P,G>& map,
	  typename fusion_property_map<P,G>::key_type key,
	  const typename fusion_property_map<P,G>::value_type& value)  {
  
    typedef fusion_property_map<P,G> map_t;
    typedef typename mpl::find<typename map_t::bundle, typename map_t::property>::type iterator;
    fusion::at_c<iterator::pos::value>(fusion::at_c<0>(map.m_graph[key])) = value;
};
  
  
template <typename P, typename G>
typename fusion_property_map<P,G>::reference at( const fusion_property_map<P,G>& map,
						 typename fusion_property_map<P,G>::key_type key)
  {
    typedef fusion_property_map<P,G> map_t;
    typedef typename mpl::find<typename map_t::bundle, typename map_t::property>::type iterator;
    fusion::at_c<iterator::pos::value>(fusion::at_c<0>(map.m_graph[key]));
  }
  
}


#endif //NS2_PROPERTY_H