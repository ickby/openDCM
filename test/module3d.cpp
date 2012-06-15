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

#include "test/Octave/debugsolver.hpp"

#include <time.h>
#include <iostream>
#include <iomanip>

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

//again, two vectors perpendicular, maybe the easiest constraints of them all
template< typename Kernel, typename Tag1, typename Tag2 >
struct test_constraint {
    typedef typename Kernel::number_type Scalar;
    typedef typename Kernel::VectorMap   Vector;

    Scalar calculate(Vector& param1,  Vector& param2)  {
        assert(false);
    };

    Scalar calculateGradientFirst(Vector& param1,  Vector& param2, Vector& dparam1) {
        assert(false);
    };

    Scalar calculateGradientSecond(Vector& param1,  Vector& param2, Vector& dparam2)  {
        assert(false);
    };

    void calculateGradientFirstComplete(Vector& param1,  Vector& param2, Vector& gradient) {
        assert(false);
    };

    void calculateGradientSecondComplete(Vector& param1,  Vector& param2, Vector& gradient) {
        assert(false);
    };
};

template< typename Kernel >
struct test_constraint<Kernel, dcm::tag::direction3D, dcm::tag::direction3D> {

    typedef typename Kernel::number_type Scalar;
    typedef typename Kernel::VectorMap   Vector;

    Scalar calculate(Vector& param1,  Vector& param2) {
        return param1.dot(param2);
    };
    Scalar calculateGradientFirst(Vector& param1, Vector& param2, Vector& dparam1) {

        return dparam1.dot(param2);
    };
    Scalar calculateGradientSecond(Vector& param1, Vector& param2, Vector& dparam2) {

        return param1.dot(dparam2);
    };
    void calculateGradientFirstComplete(Vector& param1, Vector& param2, Vector& gradient) {

        gradient(0) = param2(0);
        gradient(1) = param2(1);
        gradient(2) = param2(2);
    };
    void calculateGradientSecondComplete(Vector& param1, Vector& param2, Vector& gradient) {

        gradient(0) = param1(0);
        gradient(1) = param1(1);
        gradient(2) = param1(2);
    };
};


using namespace dcm;

BOOST_AUTO_TEST_SUITE(Module3D_test_suit);

typedef dcm::Kernel<double> Kernel;
typedef Module3D< mpl::vector<point, Eigen::Vector3d, line_t > > Module;
typedef Module3D< mpl::vector<point, Eigen::Vector3d, line_t >, std::string > ModuleID;
typedef System<Kernel, Module::type> SystemNOID;
typedef System<Kernel, ModuleID::type> SystemID;
typedef typename Module::type<SystemNOID>::Geometry3D geom;
typedef typename ModuleID::type<SystemID>::Geometry3D geomid;
typedef boost::shared_ptr<geom> geom_ptr;
typedef boost::shared_ptr<geomid> geomid_ptr;

typedef typename Module::type<SystemNOID>::Constraint3D cons;
typedef typename ModuleID::type<SystemID>::Constraint3D consid;
typedef boost::shared_ptr<cons> cons_ptr;
typedef boost::shared_ptr<consid> consid_ptr;

typedef typename SystemNOID::Cluster::vertex_iterator viter;
typedef typename Module::type<SystemNOID>::vertex_prop vertex_prop;


BOOST_AUTO_TEST_CASE(module3d_basic_solving) {

    SystemNOID sys;

    point p1,p2,p3;
    p1.push_back(7);
    p1.push_back(-0.5);
    p1.push_back(0.3);
    p2.push_back(0.2);
    p2.push_back(0.5);
    p2.push_back(-0.1);
    p3.push_back(1.2);
    p3.push_back(5.9);
    p3.push_back(0.43);


    geom_ptr g1 = sys.createGeometry3D(p1);
    geom_ptr g2 = sys.createGeometry3D(p2);
    geom_ptr g3 = sys.createGeometry3D(p3);

    //check empty solving
    sys.solve();

    //simple constraint and fire
    cons_ptr c1 = sys.createConstraint3D<test_constraint>(g1, g2);
    cons_ptr c2 = sys.createConstraint3D<test_constraint>(g2, g3);
    cons_ptr c3 = sys.createConstraint3D<test_constraint>(g3, g1);
    sys.solve();

    typename Kernel::Vector3 v1,v2,v3;
    point& rp1 = get<point>(g1);
    point& rp2 = get<point>(g2);
    point& rp3 = get<point>(g3);

    v1<<rp1[0],rp1[1],rp1[2];
    v2<<rp2[0],rp2[1],rp2[2];
    v3<<rp3[0],rp3[1],rp3[2];

    BOOST_CHECK(Kernel::isSame(v1.dot(v2),0));
    BOOST_CHECK(Kernel::isSame(v2.dot(v3),0));
    BOOST_CHECK(Kernel::isSame(v3.dot(v1),0));

}

BOOST_AUTO_TEST_CASE(module3d_cluster_solving) {

    SystemNOID sys;

    point p1,p2,p3;
    p1.push_back(7);
    p1.push_back(-0.5);
    p1.push_back(0.3);
    p2.push_back(0.2);
    p2.push_back(0.5);
    p2.push_back(-0.1);
    p3.push_back(1.2);
    p3.push_back(5.9);
    p3.push_back(0.43);

    geom_ptr g1 = sys.createGeometry3D(p1);
    geom_ptr g2 = sys.createGeometry3D(p2);
    geom_ptr g3 = sys.createGeometry3D(p3);

    //now trick a bit and move two geometries manual in a subcluster
    std::pair<typename SystemNOID::Cluster&, LocalVertex> sc = sys.m_cluster.createCluster();
    sys.m_cluster.moveToSubcluster(sys.m_cluster.getLocalVertex(g1->getProperty<vertex_prop>()).first, sc.second);
    sys.m_cluster.moveToSubcluster(sys.m_cluster.getLocalVertex(g2->getProperty<vertex_prop>()).first, sc.second);
    sc.first.setClusterProperty<changed_prop>(true);
    sc.first.setClusterProperty<type_prop>(details::cluster3D);

    //and finally add constraints
    cons_ptr c1 = sys.createConstraint3D<test_constraint>(g1, g2);
    cons_ptr c2 = sys.createConstraint3D<test_constraint>(g2, g3);
    cons_ptr c3 = sys.createConstraint3D<test_constraint>(g3, g1);

    sys.solve();

    typename Kernel::Vector3 v1,v2,v3;
    point& rp1 = get<point>(g1);
    point& rp2 = get<point>(g2);
    point& rp3 = get<point>(g3);

    v1<<rp1[0],rp1[1],rp1[2];
    v2<<rp2[0],rp2[1],rp2[2];
    v3<<rp3[0],rp3[1],rp3[2];

    BOOST_CHECK(Kernel::isSame(v1.dot(v2),0));
    BOOST_CHECK(Kernel::isSame(v2.dot(v3),0));
    BOOST_CHECK(Kernel::isSame(v3.dot(v1),0));
};

BOOST_AUTO_TEST_CASE(module3d_id) {

    SystemID sys;
    Eigen::Vector3d p1,p2;
    p1 << 7, -0.5, 0.3;
    p2 << 0.2, 0.5, -0.1;

    geomid_ptr g1 = sys.createGeometry3D(p1, "g1");
    geomid_ptr g2 = sys.createGeometry3D(p2, "g2");

    consid_ptr c1 = sys.createConstraint3D<Distance3D>("constraint", g1, g2, 5);

    BOOST_CHECK(!g1->getIdentifier().compare("g1"));
    BOOST_CHECK(!g2->getIdentifier().compare("g2"));
    BOOST_CHECK(!c1->getIdentifier().compare("constraint"));

    BOOST_CHECK(sys.hasGeometry3D("g1"));
    BOOST_CHECK(!sys.hasGeometry3D("fail"));
    BOOST_CHECK(sys.hasConstraint3D("constraint"));
    BOOST_CHECK(!sys.hasConstraint3D("fail"));

    BOOST_CHECK(sys.getGeometry3D("g1") == g1);
    BOOST_CHECK(sys.getGeometry3D("g2") == g2);
    BOOST_CHECK(sys.getConstraint3D("constraint") == c1);

}

BOOST_AUTO_TEST_CASE(module3d_distance_constraint) {

    SystemNOID sys1; //point point distance

    Eigen::Vector3d p1,p2;
    p1 << 7, -0.5, 0.3;
    p2 << 0.2, 0.5, -0.1;

    geom_ptr g1 = sys1.createGeometry3D(p1);
    geom_ptr g2 = sys1.createGeometry3D(p2);

    cons_ptr c1 = sys1.createConstraint3D<Distance3D>(g1, g2, 5);

    sys1.solve();

    Eigen::Vector3d v1,v2;
    v1 = get<Eigen::Vector3d>(g1);
    v2 = get<Eigen::Vector3d>(g2);

    BOOST_CHECK(Kernel::isSame((v1-v2).norm(), 5));

    SystemNOID sys2; //point point distance in clusters (first time to check translation)

    geom_ptr g3 = sys2.createGeometry3D(p1);
    geom_ptr g4 = sys2.createGeometry3D(p2);

    //now trick a bit and move geometries manual in a subcluster
    std::pair<typename SystemNOID::Cluster&, LocalVertex> sc1 = sys2.m_cluster.createCluster();
    std::pair<typename SystemNOID::Cluster&, LocalVertex> sc2 = sys2.m_cluster.createCluster();
    sys2.m_cluster.moveToSubcluster(sys2.m_cluster.getLocalVertex(g3->getProperty<vertex_prop>()).first, sc1.second);
    sys2.m_cluster.moveToSubcluster(sys2.m_cluster.getLocalVertex(g4->getProperty<vertex_prop>()).first, sc2.second);
    sc1.first.setClusterProperty<changed_prop>(false);//dont need to solve as it's only one point
    sc1.first.setClusterProperty<changed_prop>(false);

    cons_ptr c2 = sys2.createConstraint3D<Distance3D>(g3, g4, 5);

    sys2.solve();

    v1 = get<Eigen::Vector3d>(g3);
    v2 = get<Eigen::Vector3d>(g4);

    BOOST_CHECK(Kernel::isSame((v1-v2).norm(), 5));
}

BOOST_AUTO_TEST_CASE(module3d_parallel_constraint) {

    SystemNOID sys1; //line line parallel

    line_t l1,l2;
    l1 << 7, -0.5, 0.3, 2.3, 1.2, -0.2;
    l2 << 0.2, 0.5, -0.1, -2.1, 1.2, 0;

    geom_ptr g1 = sys1.createGeometry3D(l1);
    geom_ptr g2 = sys1.createGeometry3D(l2);

    cons_ptr c1 = sys1.createConstraint3D<Parallel3D>(g1, g2, Same);

    sys1.solve();

    line_t rl1,rl2;
    rl1 = get<line_t>(g1);
    rl2 = get<line_t>(g2);

    BOOST_CHECK(Kernel::isSame((rl1.tail<3>()-rl2.tail<3>()).norm(), 0));
}

BOOST_AUTO_TEST_CASE(module3d_angle_constraint) {

    SystemNOID sys1; //line line parallel

    line_t l1,l2;
    l1 << 7, -0.5, 0.3, 2.3, 1.2, -0.2;
    l2 << 0.2, 0.5, -0.1, -2.1, 1.2, 0;

    geom_ptr g1 = sys1.createGeometry3D(l1);
    geom_ptr g2 = sys1.createGeometry3D(l2);

    cons_ptr c1 = sys1.createConstraint3D<Angle3D>(g1, g2, 0.2);

    sys1.solve();

    line_t rl1,rl2;
    rl1 = get<line_t>(g1);
    rl2 = get<line_t>(g2);

    BOOST_CHECK(Kernel::isSame(std::acos((rl1.tail<3>().dot(rl2.tail<3>())) / (rl1.tail<3>().norm()*rl2.tail<3>().norm())) , 0.2));
}

BOOST_AUTO_TEST_SUITE_END();
