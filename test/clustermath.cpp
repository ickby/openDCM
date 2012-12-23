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

#include <boost/test/unit_test.hpp>

#include<opendcm/Core>
#include<opendcm/Module3D>

namespace dcm {

template<>
struct geometry_traits<Eigen::Vector3d> {
    typedef tag::point3D tag;
    typedef modell::XYZ modell;
    typedef orderd_roundbracket_accessor accessor;
};

}

typedef dcm::Kernel<double> Kernel;
typedef dcm::Module3D< mpl::vector< Eigen::Vector3d> > Module3D;
typedef dcm::System<Kernel, Module3D::type> System;

typedef typename Module3D::type<System>::Geometry3D Geometry3D;
typedef boost::shared_ptr<Geometry3D> Geom;

//namespace dcm {
typedef dcm::details::ClusterMath<System> cmath;
//};

BOOST_AUTO_TEST_SUITE(clustermath_suit);

BOOST_AUTO_TEST_CASE(clustermath_scaling) {

    System sys;
    cmath math;
    math.initFixMaps();

    for(int i=1; i<100; i++) {

        //add the amount of points
        for(int j=0; j<i; j++) {

            Eigen::Vector3d v = Eigen::Vector3d::Random()*100;
            Geom g(new Geometry3D(v, sys));
            //to calculate the local value which is used by scaling
            g->transformInverse(Eigen::Matrix3d::Identity(), Eigen::Vector3d::Zero(3));
            math.addGeometry(g);
        };

        //calculate the scale value for these points
        double scale = math.calculateClusterScale();

        //see if the scale value is correct
        if(i!=1) {
            for(int j=0; j<i; j++) {
                double val = (math.getGeometry()[j]->getPoint() - math.midpoint).norm();
                BOOST_REQUIRE_GE(val / scale , 0.7999);
                BOOST_REQUIRE_LE(val / scale , 1.2111);
            };
        } else BOOST_REQUIRE(scale==0);
	
	//see if we can set arbitrary bigger scale values. no hart checking as currently the alogrithm
	//is not perfect
	math.applyClusterScale(2.*scale, false);
        if(i!=1) {
            for(int j=0; j<i; j++) {
                double val = (math.getGeometry()[j]->getPoint() - math.m_shift).norm();
                BOOST_CHECK_GE(val / (2.*scale) , 0.7999);
                BOOST_CHECK_LE(val / (2.*scale) , 1.2111);
            };
        } else BOOST_REQUIRE(scale==0);

        math.clearGeometry();
        math.initFixMaps();
    }
}

BOOST_AUTO_TEST_SUITE_END();
