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
#include "opendcm/Module3D"
#include "opendcm/ModulePart"

#include "test/Octave/debugsolver.hpp"

#include <boost/test/unit_test.hpp>

typedef Eigen::Matrix<double, 6,1> plane_t;
struct place {
  Eigen::Quaterniond quat;
  Eigen::Vector3d trans;
  
  place():quat(1,2,3,4) {
    quat.normalize();
    trans<<1,2,3;
  };
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
                t.quat.w() = value; break;
            case 1:
                t.quat.x() = value; break;
            case 2:
                t.quat.y() = value; break;
            case 3:
                t.quat.z() = value; break;
            case 4:
                t.trans(0) = value; break;
            case 5:
                t.trans(1) = value; break;
            case 6:
                t.trans(2) = value; break;
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

BOOST_AUTO_TEST_SUITE(modulepart_suit);

typedef dcm::Kernel<double> Kernel;
typedef dcm::Module3D< mpl::vector< Eigen::Vector3d, plane_t> > Module3D;
typedef dcm::Module3D< mpl::vector< Eigen::Vector3d, plane_t>, std::string > Module3DID;
typedef dcm::ModulePart< mpl::vector< place > > ModulePart;
typedef dcm::ModulePart< mpl::vector< place >, std::string > ModulePartID;

typedef dcm::System<Kernel, Module3D::type, ModulePart::type> SystemNOID;
typedef dcm::System<Kernel, Module3DID::type, ModulePartID::type> SystemID;

typedef typename ModulePart::type<SystemNOID>::Part Part;
typedef typename ModulePartID::type<SystemID>::Part PartID;
typedef boost::shared_ptr<Part> Part_ptr;
typedef boost::shared_ptr<PartID> Partid_ptr;
typedef typename Module3D::type<SystemNOID>::Geometry3D Geometry3D;
typedef typename Module3DID::type<SystemID>::Geometry3D Geometry3DID;
typedef typename Module3DID::type<SystemID>::Constraint3D Constraint3DID;
typedef boost::shared_ptr<Geometry3D> Geom;
typedef boost::shared_ptr<Geometry3DID> GeomID;
typedef boost::shared_ptr<Constraint3DID> ConsID;


BOOST_AUTO_TEST_CASE(modulepart_basics) {

  Eigen::Vector3d p1,p3,p4;
  p1 << 7, -0.5, 0.3;
  p3 << -2, -1.3, -2.8;
  p4 << 0.2, -0.5, 1.2;
    
  SystemNOID sys;
  Part_ptr part1 = sys.createPart(place());
  Geom g1 = part1->addGeometry3D( p1 );
  
  Part_ptr part2 = sys.createPart(place());
  Geom g3 = part2->addGeometry3D( p3 );
  Geom g4 = part2->addGeometry3D( p4 );
  
  sys.createConstraint3D(g1,g3,dcm::distance=5);
  sys.createConstraint3D(g1,g4,dcm::distance=5);

  sys.solve();
  
  Eigen::Vector3d v1,v2,v3,v4;
  v1 = get<Eigen::Vector3d>(g1);
  v3 = get<Eigen::Vector3d>(g3);
  v4 = get<Eigen::Vector3d>(g4);

  BOOST_CHECK( Kernel::isSame((v1-v3).norm(),5.) );
  BOOST_CHECK( Kernel::isSame((v1-v4).norm(),5.) );
  
}

BOOST_AUTO_TEST_CASE(modulepart_transformations) {
  
  Eigen::Vector3d p1;
  p1 << 7, -0.5, 0.3;
  
  SystemNOID sys;
  place p;
  p.quat.x()=1; p.quat.y()=2; p.quat.z()=3; p.quat.w()=4;
  p.quat.normalize();
  Part_ptr part1 = sys.createPart(p);
  Geom g = part1->addGeometry3D(p1, dcm::Global);
  
 // BOOST_CHECK( (g->m_global-p.quat.conjugate()._transformVector(p1)).norm() < 0.0001 );  
}

BOOST_AUTO_TEST_CASE(modulepart_identifier) {
  
  Eigen::Vector3d p1,p2,p3,p4;
  p1 << 7, -0.5, 0.3;
  p2 << 0.2, 0.5, -0.1;
  
  p3 << -2, -1.3, -2.8;
  p4 << 0.2, -0.5, 1.2;
    
  SystemID sys;
  
  Partid_ptr part1 = sys.createPart(place(), "part1");
  GeomID g1 = part1->addGeometry3D( p1, "g1" );
  GeomID g2 = part1->addGeometry3D( p2, "g2" );
  
  Partid_ptr part2 = sys.createPart(place(), "part2");
  GeomID g3 = part2->addGeometry3D( p3 , "g3" );
  GeomID g4 = part2->addGeometry3D( p4 , "g4" );
  
  ConsID c1 =  sys.createConstraint3D("c1",g1,g3,dcm::distance=5);
  ConsID c2 = sys.createConstraint3D("c2",g2,g4,dcm::distance=5);
  
  BOOST_CHECK( !part1->getIdentifier().compare("part1") );
  BOOST_CHECK( !part2->getIdentifier().compare("part2") );
  BOOST_CHECK( !c1->getIdentifier().compare("c1") );
  BOOST_CHECK( !c2->getIdentifier().compare("c2") );
  
  BOOST_CHECK( part1->hasGeometry3D("g1") );
  BOOST_CHECK( !part1->hasGeometry3D("fail") );
  BOOST_CHECK( part2->hasGeometry3D("g4") );
  BOOST_CHECK( !part2->hasGeometry3D("g2") );
  
  BOOST_CHECK( sys.getPart("part1") == part1 );
  BOOST_CHECK( sys.getPart("part2") == part2 );
}
/*
BOOST_AUTO_TEST_CASE(modulepart_combined) {
  
  Eigen::Vector3d p1,p2,p3,p4, v1, v2, v3, v4;
  p1 << 7, -0.5, 0.3;
  p2 << 0.2, 0.5, -0.1;
  
  p3 << -2, -1.3, -2.8;
  p4 << 0.2, -0.5, 1.2;
    
  SystemID sys;
  
  Partid_ptr part1 = sys.createPart(place(), "part1");
  GeomID g1 = part1->addGeometry3D( p1, "g1" );
  GeomID g2 = part1->addGeometry3D( p2, "g2" );
  
  GeomID g3 = sys.createGeometry3D( p3 , "g3" );
  GeomID g4 = sys.createGeometry3D( p4 , "g4" );
  
  ConsID c1 =  sys.createConstraint3D("c1",g1,g3,dcm::distance=5);
  ConsID c2 =  sys.createConstraint3D("c2",g2,g4,dcm::distance=5);
  ConsID c3 =  sys.createConstraint3D("c3",g3,g4,dcm::distance=7);
  
  sys.solve();
  
  v1 = get<Eigen::Vector3d>(g1);
  v2 = get<Eigen::Vector3d>(g2);
  v3 = get<Eigen::Vector3d>(g3);
  v4 = get<Eigen::Vector3d>(g4);

  BOOST_CHECK(Kernel::isSame((v1-v3).norm(), 5.));
  BOOST_CHECK(Kernel::isSame((v2-v4).norm(), 5.));
  BOOST_CHECK(Kernel::isSame((v3-v4).norm(), 7.));
}

BOOST_AUTO_TEST_CASE(modulepart_fixpart) {
  
  plane_t p1,p2,p3,p4, v1, v2, v3, v4;
  p1 << 7, -0.5, 0.3, 1,0,0;
  p2 << 0.2, 0.5, -0.1, 0,1,0;  
  p3 << 0.2, -0.5, 1.2, 0,1,0;
  p4 << -2, -1.3, -2.8, 1,0,0;
  
  SystemID sys;
  
  Partid_ptr part1 = sys.createPart(place(), "part1");
  part1->fix(true);
  
  GeomID g1 = part1->addGeometry3D( p1, "g1" );
  GeomID g2 = part1->addGeometry3D( p2, "g2" );
  
  Partid_ptr part2 = sys.createPart(place(), "part2");
  GeomID g3 = part2->addGeometry3D( p3 , "g3" );
  GeomID g4 = part2->addGeometry3D( p4 , "g4" );
  
  ConsID c1 =  sys.createConstraint3D("c1",g1,g3,dcm::distance=0);
  ConsID c2 =  sys.createConstraint3D("c2",g1,g3,dcm::parallel=dcm::Same);
  ConsID c3 =  sys.createConstraint3D("c3",g2,g4,dcm::distance=0);
  ConsID c4 =  sys.createConstraint3D("c4",g2,g4,dcm::parallel=dcm::Same);
  
  sys.solve();
  
  v1 = get<plane_t>(g1);
  v2 = get<plane_t>(g2);
  v3 = get<plane_t>(g3);
  v4 = get<plane_t>(g4);
  
  BOOST_CHECK(Kernel::isSame((v1-p1).norm(),0));
  BOOST_CHECK(Kernel::isSame((v2-p2).norm(),0));
  BOOST_CHECK(Kernel::isSame((v1.tail(3)-v3.tail(3)).norm(),0));
  BOOST_CHECK(Kernel::isSame((v2.tail(3)-v4.tail(3)).norm(),0));
  BOOST_CHECK(Kernel::isSame((v1.head(3)-v3.head(3)).dot(v3.tail(3)) / v3.tail(3).norm(), 0.));
  BOOST_CHECK(Kernel::isSame((v2.head(3)-v4.head(3)).dot(v4.tail(3)) / v4.tail(3).norm(), 0.));
}

BOOST_AUTO_TEST_CASE(modulepart_idendityquaternion) {

  plane_t p1,p2;
  p1 << 7, -0.5, 0.3, 9, -2, 3;
  p2 << -2, -1.3, -2.8, 1.2, 0.3, -1.3;
  Eigen::Vector3d p3, p4;
  p3<<1,1,1;
  p4<<0,0,0;
  
  //directions should be normalized
  p1.tail<3>().normalize();
  p2.tail<3>().normalize();
 
  place p;
  p.quat.x()=0;
  p.quat.y()=0;
  p.quat.z()=0;
  p.quat.w()=1;
  p.quat.normalize();
    
  SystemNOID sys;
  Part_ptr part1 = sys.createPart(p);
  Geom g1 = part1->addGeometry3D( p1 );
  Geom g3 = part1->addGeometry3D( p3 );
  
  Part_ptr part2 = sys.createPart(p);
  Geom g2 = part2->addGeometry3D( p2 );
  Geom g4 = part2->addGeometry3D( p4 );
  
  sys.createConstraint3D(g1,g2,dcm::parallel = dcm::Same);
  sys.createConstraint3D(g3,g4,dcm::distance = 5);

  sys.solve();
  
  plane_t v1,v2;
  v1 = get<plane_t>(g1);
  v2 = get<plane_t>(g2);
  
  Eigen::Vector3d v3, v4;
  v3 = get<Eigen::Vector3d>(g3);
  v4 = get<Eigen::Vector3d>(g4);

  BOOST_CHECK( Kernel::isSame((v1.tail<3>()-v2.tail<3>()).norm(),0.) );
  BOOST_CHECK( Kernel::isSame((v3-v4).norm(),5.) );

}*/

BOOST_AUTO_TEST_SUITE_END();