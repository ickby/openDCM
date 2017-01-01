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

#include <boost/test/unit_test.hpp>
#include <boost/exception/get_error_info.hpp>

#include <Eigen/Core>

#include "opendcm/core.hpp"
#include "opendcm/module3d.hpp"
#include "derivativetest.hpp"
#include "opendcm/core/clustergraph.hpp"

typedef dcm::Eigen3Kernel<double> K;
typedef Eigen::Matrix<double, 3, 1> Vector3;

namespace dcm {
template<>
struct geometry_traits<Vector3> {
    typedef dcm::Point3 type;
    typedef dcm::modell::CartesianPoint     modell;
    typedef dcm::accessor::OrderdBracket    accessor;
};
}

typedef dcm::System<dcm::Module3D<Vector3>> System;

BOOST_AUTO_TEST_SUITE(Module3D_test_suit);
/*
BOOST_AUTO_TEST_CASE(cluster) {

    typedef dcm::numeric::Cluster3d<K>::ParameterIterator     cParIt;
    typedef dcm::numeric::Cluster3d<K>::DerivativePack        cDer;
    typedef dcm::numeric::Cluster3dGeometry<K, dcm::geometry::Point3>::DerivativePack clDer;

    try {

        dcm::numeric::LinearSystem<K> sys(10,10);
        auto cluster = std::make_shared<dcm::numeric::Cluster3d<K>>();

        cluster->init(sys);
        cluster->calculate();
                
        BOOST_CHECK(cluster->parameters().size()==6);
        BOOST_CHECK(cluster->derivatives().size()==6);
        
        //test if we change the result from the parameters
        auto p = cluster->parameters();
        for(int i=0; i<3; ++i) {
            auto initialR = cluster->rotation();
            auto initialT = cluster->translation();
            *p[i].Value = 10;
            cluster->calculate();
            BOOST_REQUIRE(!initialT.isApprox(cluster->translation()));
            BOOST_REQUIRE(initialR.isApprox(cluster->rotation()));
            BOOST_REQUIRE_EQUAL(cluster->translation()(i), 10);
        }
        for(int i=3; i<6; ++i) {
            auto initialR = cluster->rotation();
            auto initialT = cluster->translation();
            *p[i].Value = 10;
            cluster->calculate();
            BOOST_REQUIRE(initialT.isApprox(cluster->translation()));
            BOOST_REQUIRE(!initialR.isApprox(cluster->rotation()));
        }
        //reset transform
        for(auto p : cluster->parameters())
            *p.Value = 0;

        auto clGeom = std::make_shared<dcm::numeric::Cluster3dGeometry<K, dcm::geometry::Point3>>();
        clGeom->setInputEquation(cluster);
        clGeom->init(sys);
        clGeom->calculate();

        BOOST_CHECK_EQUAL(clGeom->parameters().size(),0);
        BOOST_CHECK_EQUAL(clGeom->derivatives().size(),6);

        cluster->addClusterGeometry(clGeom);
        cluster->calculate();

        for(clDer& der : clGeom->derivatives())
            BOOST_CHECK(der.second.Value != nullptr);

        //let's test the derivatives and see if we calculate them correct for 0 values
        DerivativeTest::isCorrect(std::static_pointer_cast<dcm::numeric::Equation<K, dcm::geometry::Point3<K>>>(clGeom), 
               [&]() {
                   cluster->calculate();
               }
        );
        
        //and check the derivatives for an arbitrary value
        clGeom->point() << 1,2,3;
        clGeom->point().normalize();
        clGeom->setupLocal();
        DerivativeTest::isCorrect(std::static_pointer_cast<dcm::numeric::Equation<K, dcm::geometry::Point3<K>>>(clGeom), 
               [&]() {
                   cluster->calculate();
                }
        );


    }
    catch(boost::exception& x) {
        BOOST_FAIL(*boost::get_error_info<dcm::error_message>(x));
    }
    catch(std::exception& x) {
        BOOST_FAIL("Unknown exception");
    }
};

BOOST_AUTO_TEST_CASE(geometry) {
    
    Vector3 v;
    v<<1,2,3;
    
    System s;
    std::shared_ptr<System::Geometry3D> g = s.addGeometry3D(v);
    
    BOOST_CHECK(g->holdsType());
    BOOST_CHECK(g->holdsGeometry());
    BOOST_CHECK(g->holdsGeometryType<Vector3>());
    BOOST_CHECK(g->holdsGeometryType<dcm::Point3>());
    
    Vector3 v3 = g->get<Vector3>();
    BOOST_CHECK(v3.isApprox(v));
   
    //have a look if we have the correct information in the graph
    std::shared_ptr<System::Graph> graph = std::static_pointer_cast<System::Graph>(s.getGraph());
    BOOST_CHECK(graph->vertexCount() == 1);
    BOOST_CHECK(graph->edgeCount() == 0);
    BOOST_CHECK(graph->getProperty<dcm::symbolic::GeometryProperty>(g->getVertexProperty())==g->getProperty<dcm::symbolic::GeometryProperty>());
    
    dcm::symbolic::TypeGeometry<K,dcm::geometry::Point3>* tg = static_cast<dcm::symbolic::TypeGeometry<K,dcm::geometry::Point3>*>(graph->getProperty<dcm::symbolic::GeometryProperty>(g->getVertexProperty()));
    BOOST_CHECK(tg->getPrimitveGeometry().point().isApprox(v));
    
};

BOOST_AUTO_TEST_CASE(constraint) {
      
    Vector3 v1, v2;
    v1<<1,2,3;
    v2<<4,5,6;
    
    try {
        
    System s;
    std::shared_ptr<System::Geometry3D> g1 = s.addGeometry3D(v1);
    std::shared_ptr<System::Geometry3D> g2 = s.addGeometry3D(v2);
    
    std::shared_ptr<System::Constraint3D> c1 = s.addConstraint3D(g1, g2, dcm::distance=2.);
    std::shared_ptr<System::Constraint3D> c2 = s.addConstraint3D(g1, g2, dcm::distance=3., dcm::angle=0.);
    
    //have a look if we have the correct information in the graph
    std::shared_ptr<System::Graph> graph = std::static_pointer_cast<System::Graph>(s.getGraph());
    BOOST_CHECK(graph->vertexCount() == 2);
    BOOST_CHECK(graph->edgeCount() == 1);
    std::pair<dcm::graph::LocalEdge, bool> edge = graph->edge(g1->getVertexProperty(), g2->getVertexProperty());
    BOOST_REQUIRE(edge.second);
    BOOST_CHECK(graph->getGlobalEdgeCount(edge.first) == 3);
    
    
    }
    catch(boost::exception& x) {
        BOOST_FAIL(*boost::get_error_info<dcm::error_message>(x));
    }
    catch(std::exception& x) {
        BOOST_FAIL("Unknown exception");
    }
};
*/
BOOST_AUTO_TEST_CASE(basic_solve) {
    
    Vector3 v1, v2, v3, v4;
    v1<<1,2,3;
    v2<<4,5,6;
    v3<<7,8,9;
    v4<<10,11,12;
    
    try {
        
        System s;
        s.setLoggingFilter(dcm::details::severity >= dcm::details::severity_level::iteration);
        std::shared_ptr<System::Geometry3D> g1 = s.addGeometry3D(v1);
        std::shared_ptr<System::Geometry3D> g2 = s.addGeometry3D(v2);
        std::shared_ptr<System::Geometry3D> g3 = s.addGeometry3D(v3);
        std::shared_ptr<System::Geometry3D> g4 = s.addGeometry3D(v4);
        
        std::shared_ptr<System::Constraint3D> c1 = s.addConstraint3D(g1, g2, dcm::distance=3.);
        std::shared_ptr<System::Constraint3D> c2 = s.addConstraint3D(g2, g3, dcm::distance=4.);
        std::shared_ptr<System::Constraint3D> c3 = s.addConstraint3D(g3, g4, dcm::distance=5.);
        std::shared_ptr<System::Constraint3D> c4 = s.addConstraint3D(g1, g4, dcm::distance=6.);
        
        s.solve();
    
    }
    catch(boost::exception& x) {
        std::string* message = boost::get_error_info<dcm::error_message>(x);
        if(message) 
            BOOST_FAIL(*message);
        else 
            BOOST_FAIL("Unknown boost exception");
    }
    catch(std::exception& x) {
        BOOST_FAIL("Unknown exception");
    }
}

BOOST_AUTO_TEST_SUITE_END();
