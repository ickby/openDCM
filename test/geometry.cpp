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

    BOOST_CHECK((!dcm::tag_order<test_tag1, test_tag2>::swapt::value));
    BOOST_CHECK((dcm::tag_order<test_tag2, test_tag1>::swapt::value));


    BOOST_MPL_ASSERT((boost::is_same<dcm::tag_order<test_tag1, test_tag2>::first_tag, test_tag1>));
    BOOST_MPL_ASSERT((boost::is_same<dcm::tag_order<test_tag1, test_tag2>::second_tag, test_tag2>));

    BOOST_MPL_ASSERT((boost::is_same<dcm::tag_order<test_tag2, test_tag1>::first_tag, test_tag1>));
    BOOST_MPL_ASSERT((boost::is_same<dcm::tag_order<test_tag2, test_tag1>::second_tag, test_tag2>));

}

BOOST_AUTO_TEST_CASE(geometry_transformation3d) {

    typedef dcm::Kernel<double> Kernel;

    typedef typename Kernel::Transform3D Transform;
    Transform trans3d;

    //check if initial initialisation is correct
    BOOST_CHECK(trans3d.rotation().isApprox(typename Kernel::Quaternion(1,0,0,0), 1e-10));
    BOOST_CHECK(trans3d.translation().vector().isApprox(typename Kernel::Vector3(0,0,0), 1e-10));
    BOOST_CHECK(Kernel::isSame(trans3d.scaling(), 1.));

    //check the transformations
    typename Kernel::Vector3 vec(1,2,3);
    trans3d.scale(0.5);
    vec = trans3d*vec;
    BOOST_CHECK((typename Kernel::Vector3(1,2,3)*0.5).isApprox(vec, 1e-10));

    vec << 1,2,3;
    trans3d.translate(typename Transform::Translation(1,2,3));
    trans3d.transform(vec);
    BOOST_CHECK((typename Kernel::Vector3(2,4,6)*0.5).isApprox(vec, 1e-10));

    vec << 1,2,3;
    trans3d.rotate((typename Kernel::Quaternion(1,2,3,4)).normalized());
    trans3d.transform(vec);
    typename Kernel::Vector3 res = (typename Kernel::Quaternion(1,2,3,4)).normalized()._transformVector(typename Kernel::Vector3(1,2,3));
    res += typename Kernel::Vector3(1,2,3);
    res *= 0.5;
    BOOST_CHECK(res.isApprox(vec, 1e-10));

    //check the invertion
    trans3d.invert();
    trans3d.transform(vec);
    BOOST_CHECK(vec.isApprox(typename Kernel::Vector3(1,2,3), 1e-10));

    //check successive transformations
    trans3d.setIdentity();
    trans3d *= typename Transform::Rotation(1,2,3,4);
    trans3d *= typename Transform::Translation(1,2,3);
    trans3d *= typename Transform::Scaling(2);
    Transform trans3d_2(trans3d);
    BOOST_CHECK(trans3d_2.isApprox(trans3d, 1e-10));
    BOOST_CHECK(Kernel::isSame(trans3d_2.rotation().coeffs().norm(),1));
    
    trans3d.invert();
    trans3d_2.transform(vec);
    trans3d.transform(vec);
    BOOST_CHECK(vec.isApprox(typename Kernel::Vector3(1,2,3), 1e-10));
    
    Transform trans3d_I = trans3d_2 * trans3d;
    trans3d_I.transform(vec);
    std::cout<<vec<<std::endl;
    BOOST_CHECK(vec.isApprox(typename Kernel::Vector3(1,2,3), 1e-10));

    Transform trans3d_3(typename Transform::Rotation(4,9,1,2),
            typename Transform::Translation(4,2,-6), 2);
    Transform trans3d_4(typename Transform::Rotation(1,2,4,3),
            typename Transform::Translation(-4,1,0), 3);
    Transform trans3d_5(typename Transform::Rotation(4,2,1,3),
            typename Transform::Translation(-4,-1,2), 4);
    
    vec << 1,2,3;
    trans3d_3.transform(vec);
    trans3d_4.transform(vec);
    trans3d_5.transform(vec);
    typename Kernel::Vector3 v1 = vec;

    vec << 1,2,3;
    Transform trans3d_34 = trans3d_3 * trans3d_4;
    Transform trans3d_345 = trans3d_34  * trans3d_5;
    trans3d_345.transform(vec);
    BOOST_CHECK(vec.isApprox(v1, 1e-10));
    
    vec << 1,2,3;
    trans3d_34.transform(vec);
    v1 = vec;

    trans3d_5.transform(vec);
    Transform trans3d_5I = trans3d_5.inverse();
    trans3d_5I.transform(vec);
    BOOST_CHECK(vec.isApprox(v1, 1e-10));

    vec << 1,2,3;    
    trans3d_34.transform(vec);
    v1 = vec;
    vec << 1,2,3;  
    Transform trans3d_3455I = trans3d_345 * trans3d_5I;
    trans3d_3455I.transform(vec);
    BOOST_CHECK(vec.isApprox(v1, 1e-10));
    
    vec << 1,2,3;  
    Transform trans3d_345M5 = trans3d_345 * trans3d_5.inverse();
    trans3d_345M5.transform(vec);
    BOOST_CHECK(vec.isApprox(v1, 1e-10));
    
    vec << 1,2,3;  
    trans3d_345 *= trans3d_5.inverse();
    trans3d_345.transform(vec);
    BOOST_CHECK(vec.isApprox(v1, 1e-10));
}

BOOST_AUTO_TEST_CASE(geometry_geometry3d) {
  
  
}


BOOST_AUTO_TEST_SUITE_END();
