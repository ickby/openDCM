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

#include <boost/test/unit_test.hpp>

typedef Eigen::Matrix<double, 6,1> line_t;

namespace dcm {

template<>
struct geometry_traits< Eigen::Quaterniond > {
    typedef tag::part  tag;
    typedef modell::quaternion_wxyz modell;
    typedef orderd_roundbracket_accessor accessor;
};

template<>
struct geometry_traits<Eigen::Vector3d> {
    typedef tag::point3D tag;
    typedef modell::XYZ modell;
    typedef orderd_roundbracket_accessor accessor;
};


template<>
struct geometry_traits<line_t> {
    typedef tag::line3D  tag;
    typedef modell::XYZ2 modell;
    typedef orderd_roundbracket_accessor accessor;
};

}

BOOST_AUTO_TEST_SUITE(modulepart_suit);

typedef dcm::Kernel<double> Kernel;
typedef dcm::Module3D< mpl::vector< Eigen::Vector3d, line_t> > Module3D;
typedef dcm::ModulePart< mpl::vector< Eigen::Quaterniond > > ModulePart;
typedef dcm::System<Kernel, Module3D::type, ModulePart::type> System;

typedef typename ModulePart::type<System>::Part Part;
typedef boost::shared_ptr<Part> Part_ptr;
typedef typename Module3D::type<System>::Geometry3D Geometry3D;
typedef boost::shared_ptr<Geometry3D> Geom;


BOOST_AUTO_TEST_CASE(modulepart_basics) {

  Eigen::Vector3d p1,p2,p3,p4;
  p1 << 7, -0.5, 0.3;
  p2 << 0.2, 0.5, -0.1;
  
  p3 << -2, -1.3, -2.8;
  p4 << 0.2, -0.5, 1.2;
    
  System sys;
  Part_ptr part1 = sys.createPart(Eigen::Quaterniond(1,2,5,7).normalized());
  Geom g1 = part1->addGeometry3D( p1 );
  Geom g2 = part1->addGeometry3D( p2 );
  
  Part_ptr part2 = sys.createPart(Eigen::Quaterniond(3,0.2,1,5).normalized());
  Geom g3 = part2->addGeometry3D( p3 );
  Geom g4 = part2->addGeometry3D( p4 );
  
  sys.createConstraint3D<dcm::Distance3D>(g1,g3,5);
  sys.createConstraint3D<dcm::Distance3D>(g2,g4,5);
  
  sys.solve();
  
  Eigen::Vector3d v1,v2,v3,v4;
  v1 = get<Eigen::Vector3d>(g1);
  v2 = get<Eigen::Vector3d>(g2);
  v1 = get<Eigen::Vector3d>(g3);
  v2 = get<Eigen::Vector3d>(g4);
  
  BOOST_CHECK( Kernel::isSame((v1-v3).norm(),5.) );
  BOOST_CHECK( Kernel::isSame((v2-v4).norm(),5.) );
  
}

BOOST_AUTO_TEST_SUITE_END();