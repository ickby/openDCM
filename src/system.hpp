/*
    openDCM, dimensional constraint manager
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

#ifndef NS2_SYSTEM_H
#define NS2_SYSTEM_H

#include <boost/mpl/vector.hpp>
#include <boost/mpl/map.hpp>
#include <boost/mpl/fold.hpp>
#include <boost/mpl/insert.hpp>
#include <boost/mpl/placeholders.hpp>

#include <boost/graph/adjacency_list.hpp>

#include <iostream>

#include "property.hpp"
#include "clustergraph.hpp"


namespace mpl = boost::mpl;


namespace dcm {

template<typename seq, typename state>
struct edge_fold : mpl::fold< seq, state,
            mpl::if_< is_edge_property<mpl::_2>,
            mpl::push_back<mpl::_1,mpl::_2>, mpl::_1 > > {};

template<typename seq, typename state>
struct vertex_fold : mpl::fold< seq, state,
            mpl::if_< is_vertex_property<mpl::_2>,
            mpl::push_back<mpl::_1,mpl::_2>, mpl::_1 > > {};

template<typename seq, typename state>
struct map_fold : mpl::fold< seq, state, mpl::insert<mpl::_1,mpl::_2> > {};

template<typename pair>
struct make_map : std::map<boost::shared_ptr<typename pair::first>, typename pair::second> {};

template<typename seq>
struct make_prop_maps : mpl::fold< seq, mpl::vector<>, mpl::insert<mpl::_1, make_map<mpl::_2> > > {};

struct nothing {};

template<int v>
struct empty_policy {

    template<typename T>
    struct type {
        struct inheriter {};
        typedef mpl::vector<>  	objects;
        typedef mpl::map<>	obj_propertys;
        typedef mpl::vector<>  	graph_propertys;
    };
};

template< template<class> class T1 = empty_policy<1>::type,
template<class> class T2 = empty_policy<2>::type,
template<class> class T3 = empty_policy<3>::type >
class System : 	public T1< System<T1,T2,T3> >::inheriter,
            public T2< System<T1,T2,T3> >::inheriter,
            public T3< System<T1,T2,T3> >::inheriter		{

    typedef T1< System<T1,T2,T3> > Type1;
    typedef T2< System<T1,T2,T3> > Type2;
    typedef T3< System<T1,T2,T3> > Type3;

    typedef mpl::vector<typename Type1::objects, typename Type2::objects, typename Type3::objects> objects;
    typedef mpl::map<> prop1;
    typedef typename map_fold<typename Type1::obj_propertys,
    typename map_fold<typename Type1::obj_propertys,
    typename map_fold<typename Type3::obj_propertys, prop1 >::type >::type >::type obj_propertys;

    typedef typename make_prop_maps<obj_propertys>::type obj_prop_map_vector;

    typedef mpl::vector<> edge1;
    typedef typename edge_fold< typename Type1::graph_propertys,
    typename edge_fold<typename Type2::graph_propertys,
    typename edge_fold<typename Type3::graph_propertys,edge1>::type >::type >::type edge_propertys;

    typedef mpl::vector<> vertex1;
    typedef typename vertex_fold< typename Type1::graph_propertys,
    typename vertex_fold<typename Type2::graph_propertys,
    typename vertex_fold<typename Type3::graph_propertys,vertex1>::type >::type >::type vertex_propertys;

    typedef ClusterGraph<edge_propertys, vertex_propertys, objects> Cluster;

public:
    System() {
        std::cout<<"system constructor"<<std::endl;
    };

protected:
    Cluster m_cluster;
    obj_prop_map_vector m_objMaps;


};

}

#endif //NS2_SYSTEM_H
