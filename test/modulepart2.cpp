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
#include "opendcm/modulepart.hpp"

#ifdef DCM_EXTERNAL_CORE
#include "opendcm/core/imp/system_imp.hpp"
#include "opendcm/core/imp/kernel_imp.hpp"
#endif

#ifdef DCM_EXTERNAL_3D
#include "opendcm/module3d/imp/clustermath_imp.hpp"
#include "opendcm/module3d/imp/constraint3d_imp.hpp"
#include "opendcm/module3d/imp/constraint3d_holder_imp.hpp"
#include "opendcm/module3d/imp/geometry3d_imp.hpp"
#endif

#include "test/Octave/debugsolver.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/exception/get_error_info.hpp>

typedef Eigen::Matrix<double, 6,1> plane_t;
struct place {
    Eigen::Quaterniond quat;
    Eigen::Vector3d trans;

    place():quat(1,2,3,4) {
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

namespace dcm {

template<>
struct geometry_traits< place > {
    typedef tag::part  tag;
    typedef modell::quaternion_wxyz_vec3 modell;
    typedef place_accessor accessor;
};

template<>
struct geometry_traits<Eigen::Vector3d> {
    typedef tag::point3D tag;
    typedef modell::XYZ modell;
    typedef orderd_roundbracket_accessor accessor;
};

template<>
struct geometry_traits<plane_t> {
    typedef tag::plane3D  tag;
    typedef modell::XYZ2 modell;
    typedef orderd_roundbracket_accessor accessor;
};

}

BOOST_AUTO_TEST_SUITE(modulepart2_suit);

typedef dcm::Kernel<double> Kernel;
typedef dcm::Module3D< mpl::vector2< Eigen::Vector3d, plane_t> > Module3D;
typedef dcm::Module3D< mpl::vector2< Eigen::Vector3d, plane_t>, std::string > Module3DID;
typedef dcm::ModulePart< mpl::vector1< place > > ModulePart;
typedef dcm::ModulePart< mpl::vector1< place >, std::string > ModulePartID;

typedef dcm::System<Kernel, Module3D, ModulePart> SystemNOID;
typedef dcm::System<Kernel, Module3DID, ModulePartID> SystemID;

typedef ModulePart::type<SystemNOID>::Part Part;
typedef ModulePartID::type<SystemID>::Part PartID;
typedef boost::shared_ptr<Part> Part_ptr;
typedef boost::shared_ptr<PartID> Partid_ptr;
typedef Module3D::type<SystemNOID>::Geometry3D Geometry3D;
typedef Module3DID::type<SystemID>::Geometry3D Geometry3DID;
typedef Module3DID::type<SystemID>::Constraint3D Constraint3DID;
typedef boost::shared_ptr<Geometry3D> Geom;
typedef boost::shared_ptr<Geometry3DID> GeomID;
typedef boost::shared_ptr<Constraint3DID> ConsID;


BOOST_AUTO_TEST_CASE(modulepart_cyclic) {

    try {
        Eigen::Vector3d p1,p3,p4;
        p1 << 7, -0.5, 0.3;
        p3 << -2, -1.3, -2.8;
        p4 << 0.2, -0.5, 1.2;

        SystemNOID sys;
        Part_ptr part1 = sys.createPart(place());
        Geom g1 = part1->addGeometry3D(p1);

        Part_ptr part2 = sys.createPart(place());
        Geom g3 = part2->addGeometry3D(p3);

        Part_ptr part3 = sys.createPart(place());
        Geom g4 = part3->addGeometry3D(p4);

        sys.createConstraint3D(g1,g3,dcm::distance=5.);
        sys.createConstraint3D(g3,g4,dcm::distance=5.);
        sys.createConstraint3D(g4,g1,dcm::distance=5.);

        sys.solve();

        Eigen::Vector3d v1,v3,v4;
        v1 = get<Eigen::Vector3d>(g1);
        v3 = get<Eigen::Vector3d>(g3);
        v4 = get<Eigen::Vector3d>(g4);

        BOOST_CHECK(Kernel::isSame((v1-v3).norm(),5., 1e-6));
        BOOST_CHECK(Kernel::isSame((v3-v4).norm(),5., 1e-6));
        BOOST_CHECK(Kernel::isSame((v4-v1).norm(),5., 1e-6));

    }
    catch(boost::exception& e) {
        std::cout << "Error Nr. " << *boost::get_error_info<boost::errinfo_errno>(e)
                  << ": " << *boost::get_error_info<dcm::error_message>(e)<<std::endl;
        BOOST_FAIL("Exception not expected");
    };

}

BOOST_AUTO_TEST_CASE(modulepart_rescaletransform) {

    try {
        plane_t p1,p2,p3,p4, v1, v2, v3, v4;
        p1 << 7, -0.5, 0.3, 1,0,0;
        p2 << 0.2, 0.5, -0.1, 0,1,0;
        p3 << 0.2, -0.5, 1.2, 0,1,0;
        p4 << -2, -1.3, -2.8, 1,0,0;

        SystemID sys;
        sys.setOption<dcm::solverfailure>(dcm::ApplyResults);
        sys.setOption<dcm::iterations>(0);

        place pl1, pl2;
        pl1.quat *= Kernel::Quaternion(4,3,2,1).normalized();
        Partid_ptr part1 = sys.createPart(pl1, "part1");

        GeomID g1 = part1->addGeometry3D(p1, "g1");
        GeomID g2 = part1->addGeometry3D(p2, "g2");

        Partid_ptr part2 = sys.createPart(pl2, "part2");
        GeomID g3 = part2->addGeometry3D(p3 , "g3");
        GeomID g4 = part2->addGeometry3D(p4 , "g4");

        ConsID c1 =  sys.createConstraint3D("c1",g1,g3,dcm::distance=0.);
        ConsID c2 =  sys.createConstraint3D("c2",g1,g3,dcm::orientation=dcm::equal);
        ConsID c3 =  sys.createConstraint3D("c3",g2,g4,dcm::distance=0.);
        ConsID c4 =  sys.createConstraint3D("c4",g2,g4,dcm::orientation=dcm::equal);

        sys.solve();

        place pl1t = get<place>(part1);
        place pl2t = get<place>(part2);

        BOOST_CHECK(pl1t.quat.isApprox(pl1.quat, 1e-6));
        BOOST_CHECK(pl1t.trans.isApprox(pl1.trans, 1e-6));
        BOOST_CHECK(pl2t.quat.isApprox(pl2.quat, 1e-6));
        BOOST_CHECK(pl2t.trans.isApprox(pl2.trans, 1e-6));
    }    
    catch(boost::exception& e) {
        std::cout << "Error Nr. " << *boost::get_error_info<boost::errinfo_errno>(e)
                  << ": " << *boost::get_error_info<dcm::error_message>(e)<<std::endl;
        BOOST_FAIL("Exception not expected");
    };
}

BOOST_AUTO_TEST_SUITE_END();
