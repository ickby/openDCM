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

#include "opendcm/Core"

#include <boost/test/unit_test.hpp>



struct point {};

namespace dcm {

template<>
struct geometry_traits<point> {
    typedef tag_point3D tag;
};
}

using namespace dcm;

BOOST_AUTO_TEST_SUITE(Geometry_test_suit);


BOOST_AUTO_TEST_CASE(geometry_tag_rotation) {
  

  
}


BOOST_AUTO_TEST_SUITE_END();