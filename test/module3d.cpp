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

#include <time.h>
#include <iostream>
#include <iomanip>

#include <boost/test/unit_test.hpp>

struct point : std::vector<double> {};

namespace dcm {

template<>
struct geometry_traits<point> {
    typedef tag::point3D tag;
    typedef modell::XYZ modell;
    typedef orderd_bracket_accessor accessor;
};
}

using namespace dcm;

BOOST_AUTO_TEST_SUITE(Module3D_test_suit);


BOOST_AUTO_TEST_CASE(module3d_initialising) {

    typedef dcm::Kernel<double> Kernel;
    typedef Module3D< mpl::vector<point> > Module;
    typedef System<Kernel, Module::type> System;
    System sys;
    typedef typename Module::type<System>::Geometry3D geom;
    typedef boost::shared_ptr<geom> geom_ptr;

    typedef typename Module::type<System>::Constraint3D cons;
    typedef boost::shared_ptr<cons> cons_ptr;

    geom_ptr g1 = sys.createGeometry3D(point());
    geom_ptr g2 = sys.createGeometry3D(point());

    cons_ptr c = sys.createConstraint3D<Coincident3D>(g1, g2);
    
    std::cout<<"solve system!"<<std::endl;
    sys.solve();
    
}

BOOST_AUTO_TEST_SUITE_END();
