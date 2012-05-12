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

#ifndef DCM_SYSTEM_H
#define DCM_SYSTEM_H

#include <boost/mpl/vector.hpp>
#include <boost/mpl/vector/vector0.hpp>
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

namespace details {

template<typename seq, typename state>
struct vector_fold : mpl::fold< seq, state,
            mpl::push_back<mpl::_1,mpl::_2> > {};

template<typename seq, typename state>
struct edge_fold : mpl::fold< seq, state,
            mpl::if_< is_edge_property<mpl::_2>,
            mpl::push_back<mpl::_1,mpl::_2>, mpl::_1 > > {};

template<typename seq, typename state>
struct vertex_fold : mpl::fold< seq, state,
            mpl::if_< is_vertex_property<mpl::_2>,
            mpl::push_back<mpl::_1,mpl::_2>, mpl::_1 > > {};

template<typename seq, typename state, typename obj>
struct obj_fold : mpl::fold< seq, state,
            mpl::if_< boost::is_same< details::property_kind<mpl::_2>, obj>,
            mpl::push_back<mpl::_1,mpl::_2>, mpl::_1 > > {};

template<typename objects, typename properties>
struct property_map {
    typedef typename mpl::fold<
    objects, mpl::map<>, mpl::insert< mpl::_1, mpl::pair<
    mpl::_2, details::obj_fold<properties, mpl::vector<>, mpl::_2 > > > >::type type;
};

template<typename seq, typename state>
struct map_fold : mpl::fold< seq, state, mpl::insert<mpl::_1,mpl::_2> > {};

struct nothing {};

template<int v>
struct EmptyModule {

    template<typename T>
    struct type {
        struct inheriter {};
        typedef mpl::vector<>	properties;
        typedef mpl::vector<>  objects;
    };
};

}


template< template<class> class T1 = details::EmptyModule<1>::type,
template<class> class T2 = details::EmptyModule<2>::type,
template<class> class T3 = details::EmptyModule<3>::type >
class System : 	public T1< System<T1,T2,T3> >::inheriter,
            public T2< System<T1,T2,T3> >::inheriter,
            public T3< System<T1,T2,T3> >::inheriter		{

    typedef T1< System<T1,T2,T3> > Type1;
    typedef T2< System<T1,T2,T3> > Type2;
    typedef T3< System<T1,T2,T3> > Type3;

    typedef typename details::vector_fold<typename Type3::objects,
    typename details::vector_fold<typename Type2::objects,
    typename details::vector_fold<typename Type1::objects, mpl::vector<> >::type >::type>::type objects;

    typedef typename details::vector_fold<typename Type3::properties,
    typename details::vector_fold<typename Type2::properties,
    typename details::vector_fold<typename Type1::properties, mpl::vector<> >::type >::type>::type properties;

    typedef typename details::edge_fold< properties, mpl::vector<> >::type 	edge_properties;
    typedef typename details::vertex_fold< properties, mpl::vector<> >::type 	vertex_properties;
    typedef typename details::property_map<objects, properties>::type 		object_properties;

    template<typename FT1, typename FT2, typename FT3>
    friend struct Object;

public:
    typedef ClusterGraph<edge_properties, vertex_properties, objects> Cluster;

public:
    System() {  };

    Cluster m_cluster;

};

}

#endif //DCM_SYSTEM_H

