/*
    openDCM, dimensional constraint manager
    Copyright (C) 2014  Stefan Troeger <stefantroeger@gmx.net>

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

#include <boost/test/unit_test.hpp>
#include <boost/concept_check.hpp>

#include "opendcm/core/constraint.hpp"
#include <Eigen/Core>

//two vectors perpendicular, maybe the easiest constraints of them all
struct test_constraint1 : public dcm::constraint::Constraint<int> {
    using Constraint::operator=;
    test_constraint1(const int& i) : Constraint(i) {};
    
    int& radius() {return fusion::at_c<0>(m_storage);};
};

struct test_constraint2 : public dcm::constraint::Constraint< double, char> {
    using Constraint::operator=;
    test_constraint2(const double& d, const char& c) : Constraint(d, c) {};
    
    double& radius()    {return fusion::at_c<0>(m_storage);};
    char&   direction() {return fusion::at_c<1>(m_storage);};
};

template<typename T>
void pretty(T t) {
    std::cout << __PRETTY_FUNCTION__ << std::endl;
};

BOOST_AUTO_TEST_SUITE(Constraint_test_suit);

BOOST_AUTO_TEST_CASE(module3d_multiconstraint) {
    
    test_constraint1 test1(3);
    test_constraint2 test2(0.1, 'b');
    test1 = 2;
    test2(0.2,'a');
    
    BOOST_CHECK(test1.radius() == 2);
    BOOST_CHECK(test2.radius() == 0.2);
    BOOST_CHECK(test2.direction() == 'a');
    
    test1.setDefault();
    test2.setDefault();
    
    BOOST_CHECK(test1.radius() == 3);
    BOOST_CHECK(test2.radius() == 0.1);
    BOOST_CHECK(test2.direction() == 'b');

    test1(test_constraint1(2));
    test2(test_constraint2(0.2,'a'));
     
    BOOST_CHECK(test1.radius() == 2);
    BOOST_CHECK(test2.radius() == 0.2);
    BOOST_CHECK(test2.direction() == 'a');
    
    test1 = test_constraint1(3);
    test2 = test_constraint2(0.1, 'b');
    
    BOOST_CHECK(test1.radius() == 3);
    BOOST_CHECK(test2.radius() == 0.1);
    BOOST_CHECK(test2.direction() == 'b');    
}


BOOST_AUTO_TEST_SUITE_END();
