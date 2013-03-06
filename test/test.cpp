/*
    openDCM, dimensional constraint manager
    Copyright (C) 2013  Stefan Troeger <stefantroeger@gmx.net>

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

#include "opendcm/core.hpp"
#include "opendcm/module3d.hpp"


#include <boost/test/unit_test.hpp>

struct point : std::vector<double> {};
typedef Eigen::Matrix<double, 6,1> line_t;

namespace dcm {
  
struct test_accessor {

    template<typename Scalar, int ID, typename T>
    Scalar get(T& t) {
        return 0.;
    };
    template<typename Scalar, int ID,  typename T>
    void set(Scalar value, T& t) {};
};

template<>
struct geometry_traits<point> {
    typedef tag::direction3D tag;
    typedef modell::XYZ modell;
    typedef orderd_bracket_accessor accessor;
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

template<>
struct geometry_traits<mpl::int_<1> > {
    typedef tag::point3D  tag;
    typedef modell::XYZ modell;
    typedef test_accessor accessor;
};
template<>
struct geometry_traits<mpl::int_<2> > {
    typedef tag::point3D  tag;
    typedef modell::XYZ modell;
    typedef test_accessor accessor;
};
template<>
struct geometry_traits<mpl::int_<3> > {
    typedef tag::point3D  tag;
    typedef modell::XYZ modell;
    typedef test_accessor accessor;
};
template<>
struct geometry_traits<mpl::int_<4> > {
    typedef tag::point3D  tag;
    typedef modell::XYZ modell;
    typedef test_accessor accessor;
};
template<>
struct geometry_traits<mpl::int_<6> > {
    typedef tag::point3D  tag;
    typedef modell::XYZ modell;
    typedef test_accessor accessor;
};
template<>
struct geometry_traits<mpl::int_<7> > {
    typedef tag::point3D  tag;
    typedef modell::XYZ modell;
    typedef test_accessor accessor;
};
template<>
struct geometry_traits<mpl::int_<5> > {
    typedef tag::point3D  tag;
    typedef modell::XYZ modell;
    typedef test_accessor accessor;
};
}

struct TestModule1 {

    template<typename Sys>
    struct type {
        typedef mpl::map<> signal_map;

        struct test_object1 : public dcm::Object<Sys, test_object1, signal_map > {
            test_object1(Sys& system) : dcm::Object<Sys, test_object1, signal_map >(system) { };
         };

        struct inheriter {
        };

        struct test_object1_prop {
            typedef int type;
            typedef test_object1 kind;
        };

        typedef mpl::vector1<test_object1> objects;
        typedef mpl::vector1<test_object1_prop>   properties;

        template<typename System>
        static void system_init(System& sys) {};
    };
};


typedef dcm::Kernel<double> Kernel;
typedef dcm::Module3D< mpl::vector<point, Eigen::Vector3d, line_t
, mpl::int_<1>,mpl::int_<2>,mpl::int_<3>,mpl::int_<4>,mpl::int_<5>,mpl::int_<6>,mpl::int_<7>   > > Module;
typedef dcm::System<Kernel, Module::type> System;
typedef Module::type<System>::Geometry3D geom;
typedef boost::shared_ptr<geom> geom_ptr;

typedef Module::type<System>::Constraint3D cons;
typedef boost::shared_ptr<cons> cons_ptr;

BOOST_AUTO_TEST_SUITE(test);

BOOST_AUTO_TEST_CASE(test) {

  System sys;
};

BOOST_AUTO_TEST_SUITE_END();