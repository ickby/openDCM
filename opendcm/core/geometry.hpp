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
}

namespace modell {
struct undefined {
    template<typename T, typename Res>
    void transform(T& t, Res& r) {};
};
}
struct tag_point3D  {
    typedef mpl::int_<3> parameters;
    typedef mpl::int_<1> transformations;
};
struct tag_line3D  {
    typedef mpl::int_<6> parameters;
    typedef mpl::int_<2> transformations;
};
struct tag_plane  {
    typedef mpl::int_<6> parameters;
    typedef mpl::int_<2> transformations;
};



/*
struct orderd_bracket_accessor3D {

    template<typename Kernel, typename T>
    typename Kernel::number_type get(T& t, access id) {
        return t[id];
    };
    template<typename Kernel, typename T>
    void set(typename Kernel::number_type value, T& t, access id) {
        t[id] = value;
    };
};

struct orderd_roundbracket_accessor3D {

    template<typename Scalar, int access, typename T>
    Scalar get(T& t) {
        return t(access);
    };
    template<typename Scalar, int access,  typename T>
    void set(Scalar value, T& t) {
        t(id) = value;
    };
};*/

template< typename T>
struct geometry_traits {
    typedef tag::undefined 	tag;
    typedef modell::undefined 	modell;

};

//metafunction to get the tag of a special geometry
template<typename t>
struct get_tag {
    typedef typename geometry_traits<t>::tag type;
};

}

#endif // DCM_GEOMETRY_H
