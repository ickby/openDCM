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
#include "opendcm/moduleshape3d.hpp"

namespace dcm {

template<>
struct geometry_traits<Eigen::Vector3d> {
    typedef tag::point3D tag;
    typedef modell::XYZ modell;
    typedef orderd_roundbracket_accessor accessor;
};

template<>
struct geometry_traits<Eigen::Matrix<double,6,1> > {
    typedef tag::segment3D tag;
    typedef modell::XYZ2 modell;
    typedef orderd_roundbracket_accessor accessor;
};

}

typedef dcm::Kernel<double> Kernel;
typedef dcm::Module3D< mpl::vector1<Eigen::Vector3d> > Module;
typedef dcm::System<Kernel, Module, dcm::ModuleShape3D< mpl::vector1< Eigen::Matrix<double,6,1> > > > System;
typedef Module::type<System>::Geometry3D geom;
typedef dcm::ModuleShape3D< mpl::vector1<Eigen::Matrix<double,6,1> > >::type<System>::Shape3D shape;
typedef boost::shared_ptr<geom> geom_ptr;
typedef boost::shared_ptr<shape> shape_ptr;


BOOST_AUTO_TEST_SUITE(ModuleShape3D_test_suit);

BOOST_AUTO_TEST_CASE(moduleShape3D_segment) {

    try {
        Eigen::Vector3d p1(1,2,3), p2(4,5,6);

        System sys;
        geom_ptr g1 = sys.createGeometry3D(p1);
        shape_ptr shape1 = sys.createShape3D<dcm::segment3D>(p2, g1);

        geom_ptr l = shape1->geometry(dcm::line);
        geom_ptr sp = shape1->geometry(dcm::startpoint);
        geom_ptr ep = shape1->geometry(dcm::endpoint);

        BOOST_CHECK(g1==sp || g1==ep);
        BOOST_CHECK(l->getValue().head(3).isApprox(sp->getValue(), 1e-10));
        BOOST_CHECK(l->getValue().tail(3).isApprox(ep->getValue(), 1e-10));

        Eigen::Matrix<double,6,1> s1;
        s1 << 7,8,9,10,11,12;

        shape_ptr shape2 = sys.createShape3D<dcm::segment3D>(s1);

        geom_ptr l2 = shape2->geometry(dcm::line);
        geom_ptr sp2 = shape2->geometry(dcm::startpoint);
        geom_ptr ep2 = shape2->geometry(dcm::endpoint);

	BOOST_CHECK(shape2->getValue().isApprox(s1, 1e-10));
        BOOST_CHECK(l2->getValue().isApprox(s1, 1e-10));
        BOOST_CHECK(s1.head(3).isApprox(sp2->getValue(), 1e-10));
        BOOST_CHECK(s1.tail(3).isApprox(ep2->getValue(),1e-10));

        //test the solving
        sys.createConstraint3D(sp, sp2, dcm::coincidence);
        //sys.createConstraint3D(l, l2, dcm::orientation = dcm::perpendicular );
        sys.createConstraint3D(sp, ep, dcm::distance = 5);
        sys.createConstraint3D(ep2, sp2, dcm::distance = 5);
        sys.solve();

        BOOST_CHECK(sp->getValue().isApprox(sp2->getValue(), 1e-6));
	BOOST_CHECK(shape1->getValue().isApprox(l->getValue(), 1e-10));
        BOOST_CHECK(l->getValue().head(3).isApprox(sp->getValue(), 1e-10));
        BOOST_CHECK(l->getValue().tail(3).isApprox(ep->getValue(), 1e-10));
        BOOST_CHECK(l2->getValue().isApprox(shape2->getValue(), 1e-10));
        BOOST_CHECK(shape2->getValue().head(3).isApprox(sp2->getValue(), 1e-10));
        BOOST_CHECK(shape2->getValue().tail(3).isApprox(ep2->getValue(),1e-10));

    }
    catch(boost::exception& e) {
        std::cout << "Error Nr. " << *boost::get_error_info<boost::errinfo_errno>(e)
                  << ": " << *boost::get_error_info<dcm::error_message>(e)<<std::endl;
        BOOST_FAIL("Exception not expected");
    };

};

BOOST_AUTO_TEST_SUITE_END();
