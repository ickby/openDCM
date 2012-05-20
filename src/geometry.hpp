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

struct undefined {
    typedef mpl::int_<0> parameters;
};
struct tag_point {
    typedef mpl::int_<3> parameters;
};
struct tag_line {
    typedef mpl::int_<6> parameters;
};
struct tag_plane {
    typedef mpl::int_<6> parameters;
};

template< typename T>
struct geometry_traits {
    typedef undefined tag;
};

struct Storage {
    Eigen::Vector3d point;
};

//metafunction to get the tag of a special geometry
template<typename t>
struct get_tag {
    typedef typename geometry_traits<t>::tag type;
};

}

#endif // DCM_GEOMETRY_H
