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

#include <opendcm/core.hpp>
#include <opendcm/module3d.hpp>

#include <boost/mpl/vector.hpp>

#include <boost/test/unit_test.hpp>
#include <boost/exception/get_error_info.hpp>

struct point : public Eigen::Matrix<double, 3,1> {};
struct line : public Eigen::Matrix<double, 6,1> {};
struct plane : public Eigen::Matrix<double, 6,1> {};

namespace dcm {

template<>
struct geometry_traits<point> {
    typedef tag::point3D tag;
    typedef modell::XYZ modell;
    typedef orderd_roundbracket_accessor accessor;
};

template<>
struct geometry_traits<line> {
    typedef tag::line3D  tag;
    typedef modell::XYZ2 modell;
    typedef orderd_roundbracket_accessor accessor;
};

template<>
struct geometry_traits<plane> {
    typedef tag::plane3D  tag;
    typedef modell::XYZ2 modell;
    typedef orderd_roundbracket_accessor accessor;
};

}


typedef dcm::Kernel<double> Kernel_t;
typedef dcm::Module3D< mpl::vector3<point, line, plane > > Module;
typedef dcm::System<Kernel_t, Module> System;
typedef Module::type<System>::Geometry3D geom;
typedef boost::shared_ptr<geom> geom_ptr;

typedef Module::type<System>::Constraint3D cons;
typedef boost::shared_ptr<cons> cons_ptr;

typedef System::Cluster::vertex_iterator viter;
typedef Module::type<System>::vertex_prop vertex_prop;


BOOST_AUTO_TEST_SUITE(userbug_test_suit);

BOOST_AUTO_TEST_CASE(userbug_neg_solutionspace) {

    try {
        System sys;

        point p; //point in positive x direction;
        p << 1,0,0;
        plane pl;
        pl << 0,0,0,1,0,0; //plane in x direction;

        geom_ptr g1 = sys.createGeometry3D(p);
        geom_ptr g2 = sys.createGeometry3D(pl);

        cons_ptr c = sys.createConstraint3D(g1,g2, (dcm::distance=1.) & (dcm::distance=dcm::positiv_directional));

        sys.solve();
        dcm::Distance::type<Kernel_t, dcm::tag::point3D, dcm::tag::plane3D> d;
	d.calculate(g1->getValue(), g2->getValue());
        BOOST_CHECK(d.result == 1);

        sys.removeConstraint3D(c);
        c = sys.createConstraint3D(g1,g2, (dcm::distance=1.) & (dcm::distance=dcm::negative_directional));
	sys.solve();
	
        d.calculate(g1->getValue(), g2->getValue());
        BOOST_CHECK(d.result == -1);

    }
    catch
        (boost::exception& e) {
        std::cout << "Error Nr. " << *boost::get_error_info<boost::errinfo_errno>(e)
                  << ": " << *boost::get_error_info<dcm::error_message>(e)<<std::endl;
        BOOST_FAIL("Exception not expected");
    };

};

BOOST_AUTO_TEST_SUITE_END();
