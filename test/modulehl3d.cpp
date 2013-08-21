/*
    openDCM, dimensional constraint manager
    Copyright (C) 2013  Stefan Troeger <stefantroeger@gmx.net>

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
#include <boost/exception/get_error_info.hpp>

#include "opendcm/core.hpp"
#include "opendcm/module3d.hpp"
#include "opendcm/modulehl3d.hpp"

namespace dcm {
 
template<>
struct geometry_traits<Eigen::Vector3d> {
    typedef tag::point3D tag;
    typedef modell::XYZ modell;
    typedef orderd_roundbracket_accessor accessor;
};

}

typedef dcm::Kernel<double> Kernel;
typedef dcm::Module3D< mpl::vector1<Eigen::Vector3d> > Module;
typedef dcm::System<Kernel, Module, dcm::ModuleHL3D< mpl::vector0<> > > System;
typedef Module::type<System>::Geometry3D geom;
typedef dcm::ModuleHL3D< mpl::vector0<> >::type<System>::HLGeometry3D hlgeom;
typedef boost::shared_ptr<geom> geom_ptr;
typedef boost::shared_ptr<hlgeom> hlgeom_ptr;

//test generator
struct test {

    template<typename Sys>
    struct type : public dcm::details::HLGeneratorBase<Sys> {

        //check if all needed parts are supplied
        virtual bool check() {
            
	  return true;
        };
        //initialise all relations between the geometrys, throw on error
        virtual void init() {
        };
        //get geometry3d for optional types (e.g. midpoints)
        virtual boost::shared_ptr<typename dcm::details::HLGeneratorBase<Sys>::Geometry3D> getOrCreateG3d(int type) {
	  return boost::shared_ptr<typename dcm::details::HLGeneratorBase<Sys>::Geometry3D>();
        };
        //get hlgeometry3d for optional types
        virtual boost::shared_ptr<typename dcm::details::HLGeneratorBase<Sys>::HLGeometry3D> getOrCreateHLG3d(int type) {
	  return boost::shared_ptr<typename dcm::details::HLGeneratorBase<Sys>::HLGeometry3D>();

	  
	};
    };
};


BOOST_AUTO_TEST_SUITE(ModuleHL3D_test_suit);

BOOST_AUTO_TEST_CASE(moduleHL3D_creation) {
  
  try {
      Eigen::Vector3d p1(1,2,3), p2(4,5,6);
      
  
      System sys;
      geom_ptr g1 = sys.createGeometry3D(p1);
      hlgeom_ptr hlg1 = sys.createHLGeometry3D<test>(p2, g1);
      
  }
  catch(boost::exception& e) {
    std::cout << "Error Nr. " << *boost::get_error_info<boost::errinfo_errno>(e)
      << ": " << *boost::get_error_info<dcm::error_message>(e)<<std::endl;
      BOOST_FAIL("Exception not expected");
  };
  
};

BOOST_AUTO_TEST_SUITE_END();