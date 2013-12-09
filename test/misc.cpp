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

struct point : std::vector<double> {};
typedef Eigen::Matrix<double, 6,1> line_t;

namespace dcm {

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

}
using namespace dcm;

typedef dcm::Kernel<double> Kernel_t;
typedef Module3D< mpl::vector3<point, Eigen::Vector3d, line_t > > Module;
typedef Module3D< mpl::vector3<point, Eigen::Vector3d, line_t >, std::string > ModuleID;
typedef System<Kernel_t, Module> SystemNOID;
typedef System<Kernel_t, ModuleID> SystemID;
typedef Module::type<SystemNOID>::Geometry3D geom;
typedef ModuleID::type<SystemID>::Geometry3D geomid;
typedef boost::shared_ptr<geom> geom_ptr;
typedef boost::shared_ptr<geomid> geomid_ptr;

typedef Module::type<SystemNOID>::Constraint3D cons;
typedef ModuleID::type<SystemID>::Constraint3D consid;
typedef boost::shared_ptr<cons> cons_ptr;
typedef boost::shared_ptr<consid> consid_ptr;

typedef SystemNOID::Cluster::vertex_iterator viter;
typedef Module::type<SystemNOID>::vertex_prop vertex_prop;

BOOST_AUTO_TEST_SUITE(Misc_test_suit);


BOOST_AUTO_TEST_SUITE_END();