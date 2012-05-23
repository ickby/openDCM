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


#ifndef DCM_GEOMETRY_H
#define DCM_GEOMETRY_H

#include <iostream>

#include <eigen3/Eigen/Core>

#include <boost/type_traits.hpp>
#include <boost/mpl/assert.hpp>
#include <boost/mpl/fold.hpp>
#include <boost/mpl/set.hpp>
#include <boost/mpl/less.hpp>
#include <boost/mpl/or.hpp>
#include <boost/mpl/not.hpp>
#include <boost/mpl/bool.hpp>

#include <boost/fusion/include/as_vector.hpp>
#include <boost/fusion/include/mpl.hpp>
#include <boost/concept_check.hpp>
#include <boost/graph/graph_concepts.hpp>

namespace mpl = boost::mpl;
namespace fusion = boost::fusion;

namespace dcm {

namespace tag {

struct undefined {
    typedef mpl::int_<0> parameters;
    typedef mpl::int_<0> transformations;
};

//we nee to order tags, this base make it easy for module tags
namespace weight {
struct undefined : mpl::int_<0> {};
struct point : mpl::int_<1> {};
struct line  : mpl::int_<2> {};
struct plane : mpl::int_<3> {};
}
}

namespace modell {
struct undefined {
    template<typename T, typename Res>
    void transform(T& t, Res& r) {};
};
}

struct accessor_undefined {
    template<typename Scalar, int ID, typename T>
    Scalar get(T&) {
        return ID;
    };
    template<typename Scalar, int ID, typename T>
    void set(Scalar,T&) {};
};

struct orderd_bracket_accessor {

    template<typename Scalar, int ID, typename T>
    Scalar get(T& t) {
        return t[ID];
    };
    template<typename Scalar, int ID, typename T>
    void set(Scalar value, T& t) {
        t[ID] = value;
    };
};

struct orderd_roundbracket_accessor {

    template<typename Scalar, int ID, typename T>
    Scalar get(T& t) {
        return t(ID);
    };
    template<typename Scalar, int ID,  typename T>
    void set(Scalar value, T& t) {
        t(ID) = value;
    };
};

//tag ordering
template<typename T1, typename T2>
struct tag_order {

    BOOST_MPL_ASSERT((mpl::not_< mpl::or_<
                      boost::is_same< typename T1::weight, mpl::int_<0> >,
                      boost::is_same< typename T2::weight, mpl::int_<0> >  >  >));

    typedef typename mpl::less<typename T2::weight, typename T1::weight>::type swapt;
    typedef typename mpl::if_<swapt, T2, T1>::type first_tag;
    typedef typename mpl::if_<swapt, T1, T2>::type second_tag;
};


//template<typename T1, typename T2>
//struct type_order : public tag_order< typename geometry_traits<T1>::tag, typename geometry_traits<T2>::tag > {};

template< typename T>
struct geometry_traits {
    typedef tag::undefined 	tag;
    typedef modell::undefined 	modell;
    typedef accessor_undefined  accessor;

};

}

#endif // DCM_GEOMETRY_H
