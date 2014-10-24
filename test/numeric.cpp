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

#include <boost/test/unit_test.hpp>

#include "opendcm/core/geometry.hpp"

using namespace dcm;

typedef dcm::Eigen3Kernel<double> K;

template<typename Kernel, bool MappedType = true>
struct TDirection3 : public geometry::Geometry<Kernel, MappedType,
            geometry::storage::Vector<3>> {

    using geometry::Geometry<Kernel, MappedType, geometry::storage::Vector<3>>::m_storage;
    
    auto value() -> decltype(fusion::at_c<0>(m_storage)){
        return fusion::at_c<0>(m_storage);
    };
};
 
template<typename Kernel, bool MappedType = true>
struct TLine3 : public geometry::Geometry<Kernel, MappedType,
            geometry::storage::Vector<3>, geometry::storage::Vector<3>> {

    using geometry::Geometry<Kernel, MappedType, 
        geometry::storage::Vector<3>, geometry::storage::Vector<3>>::m_storage;
    
    auto point() -> decltype(fusion::at_c<0>(m_storage)){
        return fusion::at_c<0>(m_storage);
    };
    
    auto direction() -> decltype(fusion::at_c<1>(m_storage)){
        return fusion::at_c<1>(m_storage);
    };
};
 

using namespace dcm;

template<typename T>
void pretty(T t) {
  std::cout << __PRETTY_FUNCTION__ << std::endl;  
};

BOOST_AUTO_TEST_SUITE(Numeric_test_suit);

BOOST_AUTO_TEST_CASE(geometry) {
  
    TDirection3<K>  basic;
    TDirection3<K, false> basic_vec;
    
    basic_vec.value()[0] = 5;
    
};
    
/*
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
    typedef Kernel::Transform3D Transform;
    
    Kernel k;
    Transform trans3d;

    //check if initial initialisation is correct
    BOOST_CHECK(trans3d.rotation().isApprox(Kernel::Quaternion(1,0,0,0), 1e-10));
    BOOST_CHECK(trans3d.translation().vector().isApprox(Kernel::Vector3(0,0,0), 1e-10));
    BOOST_CHECK(k.isSame(trans3d.scaling().factor(), 1.));

    //check the transformations
    Kernel::Vector3 vec(1,2,3);
    trans3d.scale(Transform::Scaling(0.5));
    vec = trans3d*vec;
    BOOST_CHECK((Kernel::Vector3(1,2,3)*0.5).isApprox(vec, 1e-10));

    vec << 1,2,3;
    trans3d.translate(Transform::Translation(1,2,3));
    trans3d.transform(vec);
    BOOST_CHECK((Kernel::Vector3(2,4,6)*0.5).isApprox(vec, 1e-10));

    vec << 1,2,3;
    trans3d.rotate((Kernel::Quaternion(1,2,3,4)).normalized());
    trans3d.transform(vec);
    Kernel::Vector3 res = (Kernel::Quaternion(1,2,3,4)).normalized()._transformVector(Kernel::Vector3(1,2,3));
    res += Kernel::Vector3(1,2,3);
    res *= 0.5;
    BOOST_CHECK(res.isApprox(vec, 1e-10));

    //check the invertion
    trans3d.invert();
    trans3d.transform(vec);
    BOOST_CHECK(vec.isApprox(Kernel::Vector3(1,2,3), 1e-10));

    //check successive transformations
    trans3d.setIdentity();
    trans3d *= Transform::Rotation(1,2,3,4);
    trans3d *= Transform::Translation(1,2,3);
    trans3d *= Transform::Scaling(2);
    Transform trans3d_2(trans3d);
    BOOST_CHECK(trans3d_2.isApprox(trans3d, 1e-10));
    BOOST_CHECK(k.isSame(trans3d_2.rotation().coeffs().norm(),1));

    trans3d.invert();
    trans3d_2.transform(vec);
    trans3d.transform(vec);
    BOOST_CHECK(vec.isApprox(Kernel::Vector3(1,2,3), 1e-10));

    Transform trans3d_I = trans3d_2 * trans3d;
    trans3d_I.transform(vec);
    BOOST_CHECK(vec.isApprox(Kernel::Vector3(1,2,3), 1e-10));

    Transform trans3d_3(Transform::Rotation(4,9,1,2),
                        Transform::Translation(4,2,-6), Transform::Scaling(2));
    Transform trans3d_4(Transform::Rotation(1,2,4,3),
                        Transform::Translation(-4,1,0), Transform::Scaling(3));
    Transform trans3d_5(Transform::Rotation(4,2,1,3),
                        Transform::Translation(-4,-1,2), Transform::Scaling(4));

    vec << 1,2,3;
    trans3d_3.transform(vec);
    trans3d_4.transform(vec);
    trans3d_5.transform(vec);
    Kernel::Vector3 v1 = vec;

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

BOOST_AUTO_TEST_CASE(geometry_geometry_link) {

    typedef dcm::Kernel<double> Kernel;

    MES mes(3,1);
    boost::shared_ptr< Geometry > g1(new Geometry()),
          g2(new Geometry());

    Eigen::Vector3d v(1,2,3);
    g1->setValue<dcm::tag::point_t>(v);
    g2->linkTo<dcm::tag::point_t>(g1, 0);

    BOOST_CHECK(g1->getValue() == g2->getValue());
}
*/

BOOST_AUTO_TEST_SUITE_END();
