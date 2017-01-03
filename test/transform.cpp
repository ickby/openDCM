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
    BOOST_CHECK(vec.isApprox(t2*Eigen::Vector3d(0,0,0)));
    
    t3.transform(vec);
    BOOST_CHECK(vec.isApprox(Eigen::Vector3d(1.,-2.,-3.)));
    BOOST_CHECK(vec.isApprox((t1*t2*t3)*Eigen::Vector3d(0.,0.,0.)));
    
    t1 *= t2*t3;
    BOOST_CHECK((t1*Eigen::Vector3d(0,0,0)).isApprox(vec));
    BOOST_CHECK((t1.inverse()*vec).isApprox(Eigen::Vector3d(0,0,0)));
    
    BOOST_CHECK(t1.invert().isApprox(t3.inverse()*t2.inverse()));
    BOOST_CHECK((t1*vec).isApprox(t3.inverse()*t2.inverse()*vec));
    BOOST_CHECK((t1*vec).isApprox(Eigen::Vector3d(0,0,0)));
    
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
    
    //check if we can assign to and can be assigned by the normal transform
    dcm::details::Transform<double,3> normal;
    t = normal;
    normal = t;   
}

BOOST_AUTO_TEST_CASE(CombineTransforms) {
    
    //the default transform types
    typedef dcm::details::Transform<double, 3>  Transform;
    Transform t1(Eigen::AngleAxisd(1*M_PI, Eigen::Vector3d::UnitX()), Eigen::Translation<double,3>(Eigen::Vector3d(1., 2., 3.)));
    Transform t2(Eigen::AngleAxisd(1.3*M_PI, Eigen::Vector3d(1.,1.,1.).normalized()), Eigen::Translation<double,3>(Eigen::Vector3d(-1., 2., -1.)));
   
    auto vec = Eigen::Vector3d(1.,1.,1.);
    t1.transform(vec);
    t2.transform(vec);
    auto t = t1*t2;
    BOOST_CHECK(vec.isApprox(t*Eigen::Vector3d(1.,1.,1.)));
    
    auto vecL = Eigen::Vector3d(1.,1.,1.);
    auto vecC = t1*vecL;
    auto vecG = t2*vecC;
    
    BOOST_CHECK(((t1*t2)*vecL).isApprox(vecG));
    BOOST_CHECK((t*vecL).isApprox(vecG));
    BOOST_CHECK(((t.inverse()*t1)*vecG).isApprox(vecC));
    BOOST_CHECK(((t2*t.inverse())*vecC).isApprox(vecL));
}


BOOST_AUTO_TEST_SUITE_END();
