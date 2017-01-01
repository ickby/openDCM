/*
    openDCM, dimensional constraint manager
    Copyright (C) 2016  Stefan Troeger <stefantroeger@gmx.net>

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

#include "opendcm/core/transformation.hpp"

#include <iostream>
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(transform_suit);

//test functions accepting all kinds of transforms 
template<typename Derived>
void testfunc(dcm::details::TransformBase<Derived>& transform) {};

BOOST_AUTO_TEST_CASE(Transform3D) {

    //the default transform types
    typedef dcm::details::Transform<double, 3>  Transform;
    typedef typename Transform::Rotation        Rotation;
    typedef typename Transform::Translation     Translation;
    
    //constructors
    Transform t1;
    Transform t2(Translation(Eigen::Vector3d(1., 2., 3.)));
    Transform t3(Rotation(Eigen::AngleAxisd(1*M_PI, Eigen::Vector3d::UnitX())));
   
    auto vec = Eigen::Vector3d(0,0,0);
    t1.transform(vec);
    BOOST_CHECK(vec.isApprox(Eigen::Vector3d(0.,0.,0.)));
    
    t2.transform(vec);
    BOOST_CHECK(vec.isApprox(Eigen::Vector3d(1.,2.,3.)));
    
    t3.transform(vec);
    BOOST_CHECK(vec.isApprox(Eigen::Vector3d(1.,-2.,-3.)));
    
    auto tn = t1*t2*t3;
    BOOST_CHECK((tn*Eigen::Vector3d(0,0,0)).isApprox(vec));
    
    t1 *= t2*t3;
    BOOST_CHECK((t1*Eigen::Vector3d(0,0,0)).isApprox(vec));
    
    testfunc(t1);
};


BOOST_AUTO_TEST_CASE(MapMatrixTransform3D) {
    
    //the default transform types
    typedef dcm::details::MapMatrixTransform<double, 3>  Transform;
    typedef Eigen::Matrix<double, 3, 3>                  Rotation;
    typedef Eigen::Matrix<double, 3, 1>                  Translation;

    //the data
    Rotation    rot = Eigen::AngleAxisd(1*M_PI, Eigen::Vector3d::UnitX()).toRotationMatrix();
    Translation trans = Eigen::Vector3d(1., 2., 3.);
    //constructors
    Transform t(rot, trans);
   
    auto vec = Eigen::Vector3d(1.,1.,1.);
    t.transform(vec);
    BOOST_CHECK(vec.isApprox(Eigen::Vector3d(2.,1.,2.)));
    
    testfunc(t);
}


BOOST_AUTO_TEST_SUITE_END();
