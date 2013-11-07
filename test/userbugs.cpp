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
#include <opendcm/modulepart.hpp>

#include <boost/mpl/vector.hpp>

#include <boost/test/unit_test.hpp>
#include <boost/exception/get_error_info.hpp>

struct point : public Eigen::Matrix<double, 3,1> {};
struct line : public Eigen::Matrix<double, 6,1> {};
struct plane : public Eigen::Matrix<double, 6,1> {};
struct cylinder : public Eigen::Matrix<double, 7,1> {};
struct place {
    Eigen::Quaterniond quat;
    Eigen::Vector3d trans;

    place() {};
    void initRandom() {
        quat = Eigen::Quaterniond(1,2,3,4);
        quat.normalize();
        trans<<1,2,3;
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
        default
                :
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

template<>
struct geometry_traits<cylinder> {
    typedef tag::cylinder3D  tag;
    typedef modell::XYZ2P modell;
    typedef orderd_roundbracket_accessor accessor;
};

template<>
struct geometry_traits< place > {
    typedef tag::part  tag;
    typedef modell::quaternion_wxyz_vec3 modell;
    typedef place_accessor accessor;
};

}


typedef dcm::Kernel<double> Kernel_t;
typedef dcm::Module3D< mpl::vector4<point, line, plane, cylinder > > Module;
typedef dcm::ModulePart< mpl::vector1<place> > ModulePart;
typedef dcm::System<Kernel_t, Module, ModulePart> System;
typedef Module::type<System>::Geometry3D geom;
typedef boost::shared_ptr<geom> geom_ptr;

typedef Module::type<System>::Constraint3D cons;
typedef boost::shared_ptr<cons> cons_ptr;

typedef ModulePart::type<System>::Part part;
typedef boost::shared_ptr<part> part_ptr;

typedef System::Cluster::vertex_iterator viter;
typedef Module::type<System>::vertex_prop vertex_prop;


BOOST_AUTO_TEST_SUITE(userbug_test_suit);
/*
BOOST_AUTO_TEST_CASE(userbug_neg_solutionspace) {

    try {
        System sys;

        point p; //point in positive x direction;
        p << 1,4,0;
        cylinder cy;
        cy << 0,0,0,1,0,0,2; //cylinder in x direction;

        geom_ptr g1 = sys.createGeometry3D(p);
        geom_ptr g2 = sys.createGeometry3D(cy);

        cons_ptr c = sys.createConstraint3D(g1,g2, (dcm::distance=1.) & (dcm::distance=dcm::positiv_directional));

        sys.solve();
        dcm::Distance::type<Kernel_t, dcm::tag::point3D, dcm::tag::cylinder3D> d;
        d.sc_value = 1;
        d.sspace = dcm::positiv_directional;
        BOOST_CHECK(d.calculate(g1->getValue(), g2->getValue())==0);
        BOOST_CHECK(d.result == 0);

        sys.removeConstraint3D(c);
        c = sys.createConstraint3D(g1,g2, (dcm::distance=1.) & (dcm::distance=dcm::negative_directional));
        sys.solve();

        d.sspace = dcm::negative_directional;
        BOOST_CHECK(d.calculate(g1->getValue(), g2->getValue())==0);
        BOOST_CHECK(d.result == -2);

        sys.removeConstraint3D(c);
        sys.removeGeometry3D(g1);
        sys.removeGeometry3D(g2);

        place pl;
        pl.initRandom();

        part_ptr p1 = sys.createPart(place());
        part_ptr p2 = sys.createPart(pl);

        g1 = p1->addGeometry3D(p);
        g2 = p2->addGeometry3D(cy);

        c = sys.createConstraint3D(g1,g2, (dcm::distance=2.) & (dcm::distance=dcm::positiv_directional));
        sys.solve();

        d.sspace = dcm::positiv_directional;
        d.sc_value=2;
        BOOST_CHECK(d.calculate(g1->getValue(), g2->getValue())==0);
        BOOST_CHECK(d.result == 0);

        sys.removeConstraint3D(c);
        c = sys.createConstraint3D(g1,g2, (dcm::distance=2.) & (dcm::distance=dcm::negative_directional));
        sys.solve();

        d.sspace = dcm::negative_directional;
        BOOST_CHECK(Kernel_t::isSame(d.calculate(g1->getValue(), g2->getValue()),0, 1e-8));
        BOOST_CHECK(Kernel_t::isSame(d.result,-4, 1e-8));

        sys.removeConstraint3D(c);
        c = sys.createConstraint3D(g1,g2, (dcm::distance=1.) & (dcm::distance=dcm::bidirectional));
        sys.solve();

        d.sc_value=1;
        d.sspace = dcm::bidirectional;
        BOOST_CHECK(Kernel_t::isSame(d.calculate(g1->getValue(), g2->getValue()),0, 1e-8));
        BOOST_CHECK(Kernel_t::isSame(d.result, -2, 1e-8));


    }
    catch
        (boost::exception& e) {
        std::cout << "Error Nr. " << *boost::get_error_info<boost::errinfo_errno>(e)
                  << ": " << *boost::get_error_info<dcm::error_message>(e)<<std::endl;
        BOOST_FAIL("Exception not expected");
    };

};*/

BOOST_AUTO_TEST_CASE(userbug_cyclic_constraints) {

    try {

        System sys;
        place p1,p2,p3;
	
	p2.quat = Eigen::Quaterniond(0.348747,-0.615122,-0.348746,0.615123);
	p2.trans << 13,-5.2512e-06,-1.01181e-05;

        //build up the system as was seen in freecad
        cylinder l1, l2, l3, l4, l5;
        l2 << -3,2.77447e-16,-4.71791e-16,-1,1.57264e-16,-9.24824e-17,1.5;
        l1 << -3,-1.5,1.75,0,0,1;
        l3 << 13,-5.2512e-06,-1.01181e-05,-1,3.28228e-07,6.32366e-07,1.5;
        l4 = l2;
        l5 << 0,3,0,0,1,0,1;

        //build up the geometry
        part_ptr pp1 = sys.createPart(p1);
        geom_ptr g1 = pp1->addGeometry3D(l1);
        geom_ptr g2 = pp1->addGeometry3D(l2);

        part_ptr pp2 = sys.createPart(p2);
        geom_ptr g3 = pp1->addGeometry3D(l3);

        part_ptr pp3 = sys.createPart(p3);
        geom_ptr g4 = pp1->addGeometry3D(l4);
        geom_ptr g5 = pp1->addGeometry3D(l5);

        //setup the constraints
	pp1->fix(true);
        sys.createConstraint3D(g1, g3, dcm::coincidence);
        sys.createConstraint3D(g2, g4, dcm::coincidence);
	sys.createConstraint3D(g3, g5, dcm::orientation=dcm::parallel);

        sys.solve();

    }
    catch
        (boost::exception& e) {
        std::cout << "Error Nr. " << *boost::get_error_info<boost::errinfo_errno>(e)
                  << ": " << *boost::get_error_info<dcm::error_message>(e)<<std::endl;
        BOOST_FAIL("Exception not expected");
    };
};

BOOST_AUTO_TEST_SUITE_END();
