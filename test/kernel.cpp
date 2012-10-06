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

#include "opendcm/Core"

#include <iostream>

#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(kernel_suit);

typedef dcm::Kernel<double> kernel;

//3 vectors constraint by 3 perpendicular constraints
struct EqnSystem : public kernel::MappedEquationSystem {

    typedef typename kernel::MappedEquationSystem Base;
    typedef typename kernel::DynStride DS;

    typename kernel::VectorMap v1, v2, v3, eqn1, eqn2, eqn3, e1_dv1, e1_dv2, e2_dv2, e2_dv3, e3_dv1, e3_dv3;

    EqnSystem() : Base(9,0,0,3), v1(NULL,0,DS(0,0)), v2(NULL,0,DS(0,0)), v3(NULL,0,DS(0,0)),
    eqn1(NULL,0,DS(0,0)), eqn2(NULL,0,DS(0,0)), eqn3(NULL,0,DS(0,0)), e1_dv1(NULL,0,DS(0,0)),
    e1_dv2(NULL,0,DS(0,0)), e2_dv2(NULL,0,DS(0,0)), e2_dv3(NULL,0,DS(0,0)),
    e3_dv3(NULL,0,DS(0,0)), e3_dv1(NULL,0,DS(0,0)) {

        int o1 = setParameterMap(dcm::Anything,3,v1);
        int o2 = setParameterMap(dcm::Anything,3,v2);
        int o3 = setParameterMap(dcm::Anything,3,v3);

        int eq1 = setResidualMap(eqn1);
        int eq2 = setResidualMap(eqn2);
        int eq3 = setResidualMap(eqn3);

        setJacobiMap(eq1,o1,3,e1_dv1);
        setJacobiMap(eq1,o2,3,e1_dv2);
        setJacobiMap(eq2,o2,3,e2_dv2);
        setJacobiMap(eq2,o3,3,e2_dv3);
        setJacobiMap(eq3,o3,3,e3_dv3);
        setJacobiMap(eq3,o1,3,e3_dv1);

    };

    virtual void recalculate() {

        eqn1 << v1.dot(v2);
        eqn2 << v2.dot(v3);
        eqn3 << v3.dot(v1);

        Base::Jacobi.setZero();
        e1_dv1 = v2;
        e1_dv2 = v1;
        e2_dv2 = v3;
        e2_dv3 = v2;
        e3_dv3 = v1;
        e3_dv1 = v3;
    };


};

struct test {
  double x,y,z;};

BOOST_AUTO_TEST_CASE(kernel_mapping) {

    typedef typename kernel::Matrix 	test_type;
    typedef typename kernel::VectorMap	test_map_type;

    test_type 	M(3,4);
    M << 1,2,3,4,5,6,7,8,9,10,11,12;

    test_map_type MM(&M(1,1), 2, typename kernel::DynStride(1,3));

    BOOST_CHECK(MM(0) == 6);
    BOOST_CHECK(MM(1) == 7);
    BOOST_CHECK(MM.size() == 2);

    test_map_type MM2(NULL,0,typename kernel::DynStride(0,0));
    BOOST_CHECK(MM2.size() == 0);

    new(&MM2) test_map_type(&M(1,1), 2, typename kernel::DynStride(1,3));

    BOOST_CHECK(MM2(0) == 6);
    BOOST_CHECK(MM2(1) == 7);
    BOOST_CHECK(MM2.size() == 2);

    typename kernel::Vector v(2);
    v << 1,2;
    MM = v;
    BOOST_CHECK(M(1,1) == 1);
    BOOST_CHECK(M(1,2) == 2);
   
    //test if fixed size works without using strides
    typename kernel::Vector3 v3;
    typename kernel::Vector3Map v3m(&v3(0)); 

    BOOST_CHECK( v3(0) == v3m(0) );
    BOOST_CHECK( v3(1) == v3m(1) );
    BOOST_CHECK( v3(2) == v3m(2) );
};

template<class T>
struct test1 {
};

BOOST_AUTO_TEST_CASE(kernel_dogleg) {

    EqnSystem s;

    s.Parameter.setRandom();
    kernel::solve(s);
    BOOST_CHECK(s.Residual.norm() < 1e-5);

}

BOOST_AUTO_TEST_SUITE_END();
