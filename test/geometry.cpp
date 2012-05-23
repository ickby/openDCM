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



struct point {
  double x,y,z;
};

namespace dcm {

//template<>
//struct geometry_traits<point> {
//    typedef tag::point3D tag;
//};
}

using namespace dcm;

BOOST_AUTO_TEST_SUITE(Geometry_test_suit);


BOOST_AUTO_TEST_CASE(geometry_accessor) {

    std::vector<double> vec;
    vec.push_back(1);
    vec.push_back(2);
    vec.push_back(3);

    BOOST_CHECK((orderd_bracket_accessor().get<double, 0>(vec) == 1));
    BOOST_CHECK((orderd_bracket_accessor().get<double, 1>(vec) == 2));
    BOOST_CHECK((orderd_bracket_accessor().get<double, 2>(vec) == 3));
    orderd_bracket_accessor().set<double, 2>(4, vec);
    BOOST_CHECK((orderd_bracket_accessor().get<double, 2>(vec) == 4));

    Eigen::Matrix<double,3,1> evec;
    evec<<1,2,3;

    BOOST_CHECK((orderd_roundbracket_accessor().get<double, 0>(evec) == 1));
    BOOST_CHECK((orderd_roundbracket_accessor().get<double, 1>(evec) == 2));
    BOOST_CHECK((orderd_roundbracket_accessor().get<double, 2>(evec) == 3));
    orderd_roundbracket_accessor().set<double, 2>(4, evec);
    BOOST_CHECK((orderd_roundbracket_accessor().get<double, 2>(evec) == 4));

}

struct test_tag1 {
  typedef dcm::tag::weight::point weight;
};
struct test_tag2 {
  typedef dcm::tag::weight::line weight;
};

BOOST_AUTO_TEST_CASE(geometry_order) {
  
  BOOST_CHECK( (!dcm::tag_order<test_tag1, test_tag2>::swapt::value) );
  BOOST_CHECK( (dcm::tag_order<test_tag2, test_tag1>::swapt::value) );
  

  BOOST_MPL_ASSERT(( boost::is_same<dcm::tag_order<test_tag1, test_tag2>::first_tag, test_tag1> ));
  BOOST_MPL_ASSERT(( boost::is_same<dcm::tag_order<test_tag1, test_tag2>::second_tag, test_tag2> ));
  
  BOOST_MPL_ASSERT(( boost::is_same<dcm::tag_order<test_tag2, test_tag1>::first_tag, test_tag1> ));
  BOOST_MPL_ASSERT(( boost::is_same<dcm::tag_order<test_tag2, test_tag1>::second_tag, test_tag2> ));

}


BOOST_AUTO_TEST_SUITE_END();
