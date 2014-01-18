/*
    openDCM, dimensional constraint manager
    Copyright (C) 2012  Stefan Troeger <stefantroeger@gmx.net>

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

#include "m3d.hpp"
#include <boost/exception/get_error_info.hpp>

test_constraint test;
comp_constraint comp_test;

using namespace dcm;

BOOST_AUTO_TEST_SUITE(Module3D_test_suit);

BOOST_AUTO_TEST_CASE(module3d_basic_solving) {

    try {

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
        cons_ptr c1 = sys.createConstraint3D(g1, g2, test);
        cons_ptr c2 = sys.createConstraint3D(g2, g3, test);
        cons_ptr c3 = sys.createConstraint3D(g3, g1, test);
        sys.solve();

        kernel::Vector3 v1,v2,v3;
        point& rp1 = get<point>(g1);
        point& rp2 = get<point>(g2);
        point& rp3 = get<point>(g3);

        v1<<rp1[0],rp1[1],rp1[2];
        v2<<rp2[0],rp2[1],rp2[2];
        v3<<rp3[0],rp3[1],rp3[2];

        BOOST_CHECK(kernel::isSame(v1.dot(v2),0, 1e-6));
        BOOST_CHECK(kernel::isSame(v2.dot(v3),0, 1e-6));
        BOOST_CHECK(kernel::isSame(v3.dot(v1),0, 1e-6));
    }
    catch(boost::exception& x) {
        BOOST_FAIL(*boost::get_error_info<error_message>(x));
    };

};

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
    std::pair<boost::shared_ptr<SystemNOID::Cluster>, LocalVertex> sc = sys.m_cluster->createCluster();
    sys.m_cluster->moveToSubcluster(sys.m_cluster->getLocalVertex(g1->getProperty<vertex_prop>()).first, sc.second);
    sys.m_cluster->moveToSubcluster(sys.m_cluster->getLocalVertex(g2->getProperty<vertex_prop>()).first, sc.second);
    sc.first->setProperty<changed_prop>(true);
    sc.first->setProperty<type_prop>(details::cluster3D);

    //and finally add constraints
    cons_ptr c1 = sys.createConstraint3D(g1, g2, test);
    cons_ptr c2 = sys.createConstraint3D(g2, g3, test);
    cons_ptr c3 = sys.createConstraint3D(g3, g1, test);

    sys.solve();

    kernel::Vector3 v1,v2,v3;
    point& rp1 = get<point>(g1);
    point& rp2 = get<point>(g2);
    point& rp3 = get<point>(g3);

    v1<<rp1[0],rp1[1],rp1[2];
    v2<<rp2[0],rp2[1],rp2[2];
    v3<<rp3[0],rp3[1],rp3[2];

    BOOST_CHECK(kernel::isSame(v1.dot(v2),0, 1e-6));
    BOOST_CHECK(kernel::isSame(v2.dot(v3),0, 1e-6));
    BOOST_CHECK(kernel::isSame(v3.dot(v1),0, 1e-6));
};

// BOOST_AUTO_TEST_CASE(module3d_multiconstraint) {
// 
//     SystemNOID sys;
// 
//     Eigen::Vector3d p1,p2,p3, p4;
//     p1 << 7, -0.5, 0.3;
//     p2 << 0.2, 0.5, -0.1;
//     p3 << -2,-1,-4;
//     p4 << 4, -3, -9;
// 
// 
//     geom_ptr g1 = sys.createGeometry3D(p1);
//     geom_ptr g2 = sys.createGeometry3D(p2);
//     geom_ptr g3 = sys.createGeometry3D(p3);
//     geom_ptr g4 = sys.createGeometry3D(p4);
// 
//     //multi constraint and fire
//     cons_ptr c1 = sys.createConstraint3D(g1, g2, test & dcm::distance(3.));
//     cons_ptr c2 = sys.createConstraint3D(g2, g3, dcm::distance(3.) & test);
//     cons_ptr c3 = sys.createConstraint3D(g3, g1, test & dcm::distance(3.) & test);
//     cons_ptr c4 = sys.createConstraint3D(g1, g4, comp_test(3.));
//     cons_ptr c5 = sys.createConstraint3D(g2, g4, comp_test(3.) & dcm::distance(3.) & comp_test(3.) & test);
//     sys.solve();
// 
//     Eigen::Vector3d& v1 = get<Eigen::Vector3d>(g1);
//     Eigen::Vector3d& v2 = get<Eigen::Vector3d>(g2);
//     Eigen::Vector3d& v3 = get<Eigen::Vector3d>(g3);
//     Eigen::Vector3d& v4 = get<Eigen::Vector3d>(g4);
// 
//     BOOST_CHECK(kernel::isSame(v1.dot(v2),0, 1e-6));
//     BOOST_CHECK(kernel::isSame((v1-v2).norm(),3, 1e-6));
//     BOOST_CHECK(kernel::isSame(v2.dot(v3),0, 1e-6));
//     BOOST_CHECK(kernel::isSame((v2-v3).norm(),3, 1e-6));
//     BOOST_CHECK(kernel::isSame(v3.dot(v1),0, 1e-6));
//     BOOST_CHECK(kernel::isSame((v1-v3).norm(),3, 1e-6));
//     BOOST_CHECK(kernel::isSame(v4.dot(v1),0, 1e-6));
//     BOOST_CHECK(kernel::isSame((v1-v4).norm(),3, 1e-6));
//     BOOST_CHECK(kernel::isSame(v4.dot(v2),0, 1e-6));
//     BOOST_CHECK(kernel::isSame((v4-v2).norm(),3, 1e-6));
// }

BOOST_AUTO_TEST_CASE(module3d_id) {

    SystemID sys;
    Eigen::Vector3d p1,p2;
    p1 << 7, -0.5, 0.3;
    p2 << 0.2, 0.5, -0.1;

    geomid_ptr g1 = sys.createGeometry3D(p1, "g1");
    geomid_ptr g2 = sys.createGeometry3D(p2, "g2");

    consid_ptr c1 = sys.createConstraint3D("constraint", g1, g2, test);

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

    BOOST_REQUIRE_THROW(sys.createConstraint3D("constraint2", g1, g2, orientation), constraint_error);
}

BOOST_AUTO_TEST_CASE(module3d_cloning) {

    SystemID sys;

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


    geomid_ptr g1 = sys.createGeometry3D(p1, "g1");
    geomid_ptr g2 = sys.createGeometry3D(p2, "g2");
    geomid_ptr g3 = sys.createGeometry3D(p3, "g3");

    //simple constraint
    consid_ptr c1 = sys.createConstraint3D("c1", g1, g2, test);
    consid_ptr c2 = sys.createConstraint3D("c2", g2, g3, test);
    consid_ptr c3 = sys.createConstraint3D("c3", g3, g1, test);

    //clone and change initial system
    SystemID* clone = sys.clone();
    p1.clear();
    p1.push_back(1);
    p1.push_back(2);
    p1.push_back(3);
    g1->set(p1);

    //check if the cloned system was affekted
    geomid_ptr cg1 = clone->getGeometry3D("g1");
    point& cp1 = get<point>(cg1);
    BOOST_CHECK(cp1[0] == 7);
    BOOST_CHECK(cp1[1] == -0.5);
    BOOST_CHECK(cp1[2] == 0.3);

    //solve and see what happens
    clone->solve();

    kernel::Vector3 v1,v2,v3;
    point& rp1 = get<point>(clone->getGeometry3D("g1"));
    point& rp2 = get<point>(clone->getGeometry3D("g2"));
    point& rp3 = get<point>(clone->getGeometry3D("g3"));

    v1<<rp1[0],rp1[1],rp1[2];
    v2<<rp2[0],rp2[1],rp2[2];
    v3<<rp3[0],rp3[1],rp3[2];

    //check if the system was solved correctly
    BOOST_CHECK(kernel::isSame(v1.dot(v2),0, 1e-6));
    BOOST_CHECK(kernel::isSame(v2.dot(v3),0, 1e-6));
    BOOST_CHECK(kernel::isSame(v3.dot(v1),0, 1e-6));

    //check if the original system is unchanged
    BOOST_CHECK(p1[0] == get<point>(sys.getGeometry3D("g1"))[0]);
    BOOST_CHECK(p1[1] == get<point>(sys.getGeometry3D("g1"))[1]);
    BOOST_CHECK(p1[2] == get<point>(sys.getGeometry3D("g1"))[2]);
    BOOST_CHECK(p2[0] == get<point>(sys.getGeometry3D("g2"))[0]);
    BOOST_CHECK(p2[1] == get<point>(sys.getGeometry3D("g2"))[1]);
    BOOST_CHECK(p2[2] == get<point>(sys.getGeometry3D("g2"))[2]);
    BOOST_CHECK(p3[0] == get<point>(sys.getGeometry3D("g3"))[0]);
    BOOST_CHECK(p3[1] == get<point>(sys.getGeometry3D("g3"))[1]);
    BOOST_CHECK(p3[2] == get<point>(sys.getGeometry3D("g3"))[2]);

    delete clone;

};

BOOST_AUTO_TEST_SUITE_END();
