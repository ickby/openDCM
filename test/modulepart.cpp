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

#include "opendcm/core.hpp"
#include "opendcm/module3d.hpp"
#include "opendcm/modulepart3d.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/exception/get_error_info.hpp>

struct Place {
    Eigen::Quaterniond quat;
    Eigen::Vector3d trans;

    Place():quat(1,2,3,4) {
        quat.normalize();
        trans<<1,2,3;
    };

    Eigen::Vector3d transformPoint(const Eigen::Vector3d& v) {
        return quat._transformVector(v) + trans;
    };

public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

struct place_accessor {

    template<typename Scalar, int ID, typename T>
    Scalar get(T& t) {
        switch(ID) {
        case 0:
            return t.quat.w();
        case 1:
            return t.quat.x();
        case 2:
            return t.quat.y();
        case 3:
            return t.quat.z();
        case 4:
            return t.trans(0);
        case 5:
            return t.trans(1);
        case 6:
            return t.trans(2);
        default:
            return 0;
        };
    };
    template<typename Scalar, int ID, typename T>
    void set(Scalar value, T& t) {
        switch(ID) {
        case 0:
            t.quat.w() = value;
            break;
        case 1:
            t.quat.x() = value;
            break;
        case 2:
            t.quat.y() = value;
            break;
        case 3:
            t.quat.z() = value;
            break;
        case 4:
            t.trans(0) = value;
            break;
        case 5:
            t.trans(1) = value;
            break;
        case 6:
            t.trans(2) = value;
            break;
        };
    };
    template<typename T>
    void finalize(T& t) {};
};

typedef dcm::Eigen3Kernel<double> K;
typedef Eigen::Matrix<double, 3, 1> Vector3;

namespace dcm {
template<>
struct geometry_traits<Vector3> {
    typedef dcm::Point3                  type;
    typedef dcm::modell::CartesianPoint  modell;
    typedef dcm::accessor::OrderdBracket accessor;
};

template<>
struct geometry_traits< Place > {
    typedef dcm::Part3                   type;
    typedef modell::quaternion_wxyz_vec3 modell;
    typedef place_accessor               accessor;
};
}

BOOST_AUTO_TEST_SUITE(ModulePart3D_test_suit);

typedef dcm::System<dcm::Module3D<Vector3>, dcm::ModulePart3D<Place>> System;


BOOST_AUTO_TEST_CASE(modulepart_basics) {

    try {
        Eigen::Vector3d p1,p3,p4;
        p1 << 7, -0.5, 0.3;
        p3 << -2, -1.3, -2.8;
        p4 << 0.2, -0.5, 1.2;

        System sys;
        std::shared_ptr<System::Part3D> part1 = sys.addPart3D(Place());
        std::shared_ptr<System::Geometry3D> g1 = part1->addGeometry3D(p1);

        std::shared_ptr<System::Part3D> part2 = sys.addPart3D(Place());
        std::shared_ptr<System::Geometry3D> g3 = part2->addGeometry3D(p3);
        std::shared_ptr<System::Geometry3D> g4 = part2->addGeometry3D(p4);

        sys.addConstraint3D(g1, g3, dcm::distance=5.);
        sys.addConstraint3D(g1, g4, dcm::distance=5.);

        sys.solve();

        Eigen::Vector3d v1,v2,v3,v4;
        v1 = g1->get<Eigen::Vector3d>();
        v3 = g3->get<Eigen::Vector3d>();
        v4 = g4->get<Eigen::Vector3d>();

        BOOST_CHECK_EQUAL((v1-v3).norm(),5.);
        BOOST_CHECK_EQUAL((v1-v4).norm(),5.);

    }
    catch(boost::exception& e) {
        std::cout << "Error Nr. " << *boost::get_error_info<boost::errinfo_errno>(e)
                  << ": " << *boost::get_error_info<dcm::error_message>(e)<<std::endl;
        BOOST_FAIL("Exception not expected");
    };
}


BOOST_AUTO_TEST_SUITE_END();
