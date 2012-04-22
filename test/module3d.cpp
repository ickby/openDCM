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

#include "system.hpp"
#include "module3d.hpp"

//#define BOOST_TEST_MODULE Geometry3DModule
#include <boost/test/unit_test.hpp>

using namespace dcm;

BOOST_AUTO_TEST_SUITE( Module3D_test_suit);

BOOST_AUTO_TEST_CASE(initialising) {
  
  typedef System<Module3D::type> System;
  System sys;
  typedef typename Module3D::type<System>::Geometry3D geom;
  typedef boost::shared_ptr<geom> geom_ptr;
 
  geom_ptr p = sys.createGeometry3D();
  
}

BOOST_AUTO_TEST_SUITE_END();