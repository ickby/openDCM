/*
    openDCM, dimensional constraint manager
    Copyright (C) 2012  Stefan Troeger <stefantroeger@gmx.net>

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

#include "opendcm/core/kernel.hpp"

#include <iostream>

#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(kernel_suit);

typedef dcm::Kernel<double> kernel;

//3 vectors constraint by 3 perpendicular constraints
struct EqnSystem : public kernel::MappedEquationSystem {

    typedef kernel::MappedEquationSystem Base;
    typedef kernel::DynStride DS;

    kernel::VectorMap v1, v2, v3, eqn1, eqn2, eqn3, e1_dv1, e1_dv2, e2_dv2, e2_dv3, e3_dv1, e3_dv3;

    EqnSystem() : Base(9,3), v1(NULL,0,DS(0,0)), v2(NULL,0,DS(0,0)), v3(NULL,0,DS(0,0)),
    eqn1(NULL,0,DS(0,0)), eqn2(NULL,0,DS(0,0)), eqn3(NULL,0,DS(0,0)), e1_dv1(NULL,0,DS(0,0)),
    e1_dv2(NULL,0,DS(0,0)), e2_dv2(NULL,0,DS(0,0)), e2_dv3(NULL,0,DS(0,0)),
    e3_dv3(NULL,0,DS(0,0)), e3_dv1(NULL,0,DS(0,0)) {

        int o1 = setParameterMap(3,v1);
        int o2 = setParameterMap(3,v2);
        int o3 = setParameterMap(3,v3);

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

    virtual void removeLocalGradientZeros() {};

};

struct test {
  double x,y,z;};

BOOST_AUTO_TEST_CASE(kernel_mapping) {

    typedef kernel::Matrix 	test_type;
    typedef kernel::VectorMap	test_map_type;

    test_type 	M(3,4);
    M << 1,2,3,4,5,6,7,8,9,10,11,12;

    test_map_type MM(&M(1,1), 2, kernel::DynStride(1,3));

    BOOST_CHECK(MM(0) == 6);
    BOOST_CHECK(MM(1) == 7);
    BOOST_CHECK(MM.size() == 2);

    test_map_type MM2(NULL,0, kernel::DynStride(0,0));
    BOOST_CHECK(MM2.size() == 0);

    new(&MM2) test_map_type(&M(1,1), 2, kernel::DynStride(1,3));

    BOOST_CHECK(MM2(0) == 6);
    BOOST_CHECK(MM2(1) == 7);
    BOOST_CHECK(MM2.size() == 2);

    kernel::Vector v(2);
    v << 1,2;
    MM = v;
    BOOST_CHECK(M(1,1) == 1);
    BOOST_CHECK(M(1,2) == 2);
   
    //test if fixed size works without using strides
    kernel::Vector3 v3(1,2,3);
    kernel::Vector3Map v3m(&v3(0)); 

    BOOST_CHECK( v3(0) == v3m(0) );
    BOOST_CHECK( v3(1) == v3m(1) );
    BOOST_CHECK( v3(2) == v3m(2) );
};

BOOST_AUTO_TEST_CASE(kernel_multimap) {
  
    // test all constructors
    // ********************* 
    kernel::Vector3 vec(1,2,3);    
    dcm::details::MultiMap<kernel::Vector3> map(vec); 
    BOOST_REQUIRE(map.rows() == 3);
    BOOST_REQUIRE(map.cols() == 1);
    BOOST_CHECK(map==vec);
    BOOST_CHECK(map+vec == 2*map);    
    
    kernel::Vector vec2(4,1);
    vec2 << 1,2,3,4;
    dcm::details::MultiMap<kernel::Vector> map2(vec2);    
    BOOST_REQUIRE(map2.rows() == 4);
    BOOST_REQUIRE(map2.cols() == 1);
    BOOST_CHECK(map2==vec2);
    BOOST_CHECK(map2+vec2 == 2*map2);
    
    kernel::Vector vec_comp(7);
    vec_comp << 1,2,3,4,1,2,3;
    map2.extend(vec.data(), 3, kernel::DynStride(3,1));
    BOOST_REQUIRE(map2.rows() == 7);
    BOOST_REQUIRE(map2.cols() == 1);
    BOOST_CHECK(map2==vec_comp);
    BOOST_CHECK(map2+vec_comp == 2*map2);
    
    kernel::Vector vec_comp2(10);
    vec_comp2 << 1,2,3,4,1,2,3,1,2,3;
    map2.extend(vec.data(), 3);
    BOOST_REQUIRE(map2.rows() == 10);
    BOOST_REQUIRE(map2.cols() == 1);
    BOOST_CHECK(map2==vec_comp2);
    BOOST_CHECK(map2+vec_comp2 == 2*map2);
    
    kernel::Vector vec_comp3(13);
    vec_comp3 << 1,2,3,4,1,2,3,1,2,3,1,2,3;
    BOOST_REQUIRE(map2.rows() == 13);
    BOOST_REQUIRE(map2.cols() == 1);
    BOOST_CHECK(map2==vec_comp3);
    BOOST_CHECK(map2+vec_comp3 == 2*map2);
    
    kernel::Matrix3 mat;
    mat << 1,2,3,4,5,6,7,8,9;
    dcm::details::MultiMap<kernel::Matrix3> map3(mat);    
    BOOST_REQUIRE(map3.rows() == 3);
    BOOST_REQUIRE(map3.cols() == 3);
    BOOST_CHECK(map3==mat);
    BOOST_CHECK(map3+mat == 2*mat);
    
    kernel::Matrix mat2(5,3);
    mat2.setRandom();
    dcm::details::MultiMap<kernel::Matrix> map4(mat2);    
    BOOST_REQUIRE(map4.rows() == 5);
    BOOST_REQUIRE(map4.cols() == 3);
    BOOST_CHECK(map4==mat2);
    BOOST_CHECK(map4+mat2 == 2*map4);
    
    kernel::Matrix39 mat3;
    mat3.setRandom();
    dcm::details::MultiMap<kernel::Matrix39> map5(mat3.data());    
    BOOST_REQUIRE(map5.rows() == 3);
    BOOST_REQUIRE(map5.cols() == 9);
    BOOST_CHECK(map5==mat3);
    BOOST_CHECK(map5+mat3 == 2*map5);
    
    dcm::details::MultiMap<kernel::Vector> map6(vec2.data(), vec2.rows());    
    BOOST_REQUIRE(map6.rows() == 4);
    BOOST_REQUIRE(map6.cols() == 1);
    BOOST_CHECK(map6==vec2);
    BOOST_CHECK(map6+vec2 == 2*map6);
    
    dcm::details::MultiMap<kernel::Matrix> map7(mat3.data(), mat3.rows(), mat3.cols());    
    BOOST_REQUIRE(map7.rows() == 3);
    BOOST_REQUIRE(map7.cols() == 9);
    BOOST_CHECK(map7==mat3);
    BOOST_CHECK(map7+mat3 == 2*map7);
};

template<class T>
struct test1 {
};

BOOST_AUTO_TEST_CASE(kernel_dogleg) {

    EqnSystem s;
    kernel k;

    s.Parameter.setRandom();
    k.solve(s);
    BOOST_CHECK(s.Residual.norm() < 1e-5);

}

BOOST_AUTO_TEST_CASE(kernel_transformation) {

    dcm::detail::Transform<double, 3> T;
    std::cout<<T<<std::endl;

}

BOOST_AUTO_TEST_SUITE_END();
