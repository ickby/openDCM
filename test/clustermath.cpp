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
/*
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
                BOOST_CHECK_GE(val / scale , 0.7999);
                BOOST_CHECK_LE(val / scale , 1.2111);
            };
        } else BOOST_REQUIRE(scale==0);

        //see if we can set arbitrary bigger scale values. no hart checking as currently the alogrithm
        //is not perfect
        math.applyClusterScale(2.*scale, false);
        if(i!=1) {
            for(int j=0; j<i; j++) {
                double val = (math.getGeometry()[j]->getPoint() - math.m_shift).norm();
                BOOST_CHECK_GE(val / scale , 2.*0.7999);
                BOOST_CHECK_LE(val / scale , 2.*1.2111);
            };
        } else BOOST_REQUIRE(scale==0);

        math.clearGeometry();
        math.initFixMaps();
    }
}*/

BOOST_AUTO_TEST_CASE(clustermath_identityhandling) {

    System sys;
    cmath math;
    
    typename Kernel::Quaternion Q(1,2,3,4);
    Q.normalize();
    
    //need to init a few things
    typename Kernel::Vector3 vec(0,0,0);
    math.initFixMaps();
    new(&math.m_normQ) typename Kernel::Vector3Map(&vec(0));
    math.m_original_translation.setZero();
    math.m_quaternion = Q;

    //add two points to the clustermath
    Eigen::Vector3d p1 = Eigen::Vector3d::Random()*100;
    Eigen::Vector3d p2 = Eigen::Vector3d::Random()*100;
    Geom g1(new Geometry3D(p1, sys));
    Geom g2(new Geometry3D(p2, sys));
    //the stuff that is normaly done by map downstream geometry
    math.setRotationMap(g1->getRotationMap(), g1->getDiffRotationMap());
    math.setTranslationMap(g1->getTranslationMap());
    math.setShiftMap(g1->getShiftMap());
    g1->m_offset = math.getParameterOffset();
    g1->transformInverse(Q.conjugate().toRotationMatrix(), -math.m_translation);
    g1->setClusterMode(true, false);
    math.setRotationMap(g2->getRotationMap(), g2->getDiffRotationMap());
    math.setTranslationMap(g2->getTranslationMap());
    math.setShiftMap(g2->getShiftMap());
    g2->m_offset = math.getParameterOffset();
    g2->transformInverse(Q.conjugate().toRotationMatrix(), -math.m_translation);
    g2->setClusterMode(true, false);
    math.addGeometry(g1);
    math.addGeometry(g2);
 
    //check if we have the local values right
    BOOST_CHECK(Kernel::isSame((g1->m_toplocal-Q.conjugate()._transformVector(p1)).norm(),0.));
    BOOST_CHECK(Kernel::isSame((g2->m_toplocal-Q.conjugate()._transformVector(p2)).norm(),0.));
   
    math.resetClusterRotation(Q);    

    //check if the new toplocal is changed to the new rotation and that Q is adjusted
    BOOST_CHECK(Kernel::isSame((g1->m_toplocal-Q.conjugate()._transformVector(p1)).norm(),0.));
    BOOST_CHECK(Kernel::isSame((g2->m_toplocal-Q.conjugate()._transformVector(p2)).norm(),0.));

    BOOST_CHECK(Kernel::isSame((Q._transformVector(g1->m_toplocal)-p1).norm(),0.));
    BOOST_CHECK(Kernel::isSame((Q._transformVector(g2->m_toplocal)-p2).norm(),0.));
    BOOST_CHECK(!Kernel::isSame(math.m_normQ.norm(), 0));
    
    //see if the downstream processing works
    math.m_rotation = Q.toRotationMatrix();
    g1->recalculate(1.);
    g2->recalculate(1.);
    BOOST_CHECK(Kernel::isSame((g1->m_rotated-p1).norm(),0.));
    BOOST_CHECK(Kernel::isSame((g2->m_rotated-p2).norm(),0.));
   
    //this should redo the rotation to original
    math.resetClusterRotation(Q);
    math.m_rotation = Q.toRotationMatrix();
    g1->recalculate(1.);
    g2->recalculate(1.);
    
    //let's check if rotations still hold
    BOOST_CHECK(Kernel::isSame((g1->m_toplocal-Q.conjugate()._transformVector(p1)).norm(),0.));
    BOOST_CHECK(Kernel::isSame((g2->m_toplocal-Q.conjugate()._transformVector(p2)).norm(),0.));
    
    //see if we have the same Quaternion as at the begining
    typename Kernel::Quaternion Qinit(1,2,3,4);
    Qinit.normalize();
    BOOST_CHECK(Qinit.isApprox(Q, 0.00001));
   
    //see if it works with shifting and scaling
    typename Kernel::number_type s = math.calculateClusterScale();
    math.applyClusterScale(s, false);
    math.resetClusterRotation(Q);
	
    math.recalculate();
    g1->recalculate(1./s);
    g2->recalculate(1./s);
    BOOST_CHECK(Kernel::isSame((g1->m_rotated*s-p1).norm(),0.));
    BOOST_CHECK(Kernel::isSame((g2->m_rotated*s-p2).norm(),0.));
    
    math.finishCalculation();
    g1->recalculate(1.);
    g2->recalculate(1.);
    BOOST_CHECK(Kernel::isSame((g1->m_rotated-p1).norm(),0.));
    BOOST_CHECK(Kernel::isSame((g2->m_rotated-p2).norm(),0.));
    BOOST_CHECK(math.m_quaternion.isApprox(Qinit, 0.0001));
    
    
    //see if it works with initial translation
    Q.setIdentity();
    math.m_quaternion.setIdentity();
    math.m_normQ.setZero();
    math.initFixMaps();
    math.m_translation<<1,2,3;
    g1->transformInverse(Eigen::Matrix3d::Identity(), -math.m_translation);
    g2->transformInverse(Eigen::Matrix3d::Identity(), -math.m_translation);
    
    s = math.calculateClusterScale();
    math.applyClusterScale(s, false);
    math.recalculate();
    g1->recalculate(1./s);
    g2->recalculate(1./s);
    BOOST_CHECK(Kernel::isSame((g1->m_rotated*s-p1).norm(),0.));
    BOOST_CHECK(Kernel::isSame((g2->m_rotated*s-p2).norm(),0.));
    
    math.finishCalculation();
    g1->recalculate(1.);
    g2->recalculate(1.);
    BOOST_CHECK(Kernel::isSame((g1->m_rotated-p1).norm(),0.));
    BOOST_CHECK(Kernel::isSame((g2->m_rotated-p2).norm(),0.));
}

BOOST_AUTO_TEST_SUITE_END();
