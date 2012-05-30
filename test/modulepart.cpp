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


BOOST_AUTO_TEST_CASE(modulepart_basics) {

  System sys;
  Part p(Eigen::Quaterniond(), sys, sys.m_cluster);
  
}

BOOST_AUTO_TEST_SUITE_END();