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

#include "parser.hpp"


BOOST_AUTO_TEST_SUITE(parser_suit);

BOOST_AUTO_TEST_CASE(parser_graph) {

    System sys; 

    fusion::vector<dcm::LocalVertex, dcm::GlobalVertex> res1 = sys.m_cluster.addVertex();
    fusion::vector<dcm::LocalVertex, dcm::GlobalVertex> res2 = sys.m_cluster.addVertex();
    fusion::vector<dcm::LocalVertex, dcm::GlobalVertex> res3 = sys.m_cluster.addVertex();
    dcm::LocalEdge e1 = fusion::at_c<0>(sys.m_cluster.addEdge(fusion::at_c<1>(res1),fusion::at_c<1>(res2)));
    dcm::LocalEdge e2 = fusion::at_c<0>(sys.m_cluster.addEdge(fusion::at_c<1>(res2),fusion::at_c<1>(res3)));
    dcm::LocalEdge e3 = fusion::at_c<0>(sys.m_cluster.addEdge(fusion::at_c<1>(res1),fusion::at_c<1>(res3)));

    boost::shared_ptr<TestModule1::type<System>::test_object1> ptr(new TestModule1::type<System>::test_object1(sys));
    sys.m_cluster.setObject(fusion::at_c<1>(res1), ptr);  
    
    sys.m_cluster.setProperty<TestModule1::type<System>::test_edge1_prop>(e1, 1);
    sys.m_cluster.setProperty<TestModule1::type<System>::test_edge1_prop>(e2, 2);
    sys.m_cluster.setProperty<TestModule1::type<System>::test_edge1_prop>(e3, 3);  
   
    //subcluster
    System::Cluster& scl1 = sys.m_cluster.createCluster().first;
    System::Cluster& scl2 = sys.m_cluster.createCluster().first;
    
    fusion::vector<dcm::LocalVertex, dcm::GlobalVertex> res4 = scl1.addVertex();
    fusion::vector<dcm::LocalVertex, dcm::GlobalVertex> res5 = scl1.addVertex();
    fusion::vector<dcm::LocalVertex, dcm::GlobalVertex> res6 = scl2.addVertex();
    dcm::LocalEdge e4 = fusion::at_c<0>(sys.m_cluster.addEdge(fusion::at_c<1>(res4),fusion::at_c<1>(res5)));
    dcm::LocalEdge e5 = fusion::at_c<0>(sys.m_cluster.addEdge(fusion::at_c<1>(res5),fusion::at_c<1>(res6)));
    dcm::LocalEdge e6 = fusion::at_c<0>(sys.m_cluster.addEdge(fusion::at_c<1>(res1),fusion::at_c<1>(res4)));
    
    std::stringstream s;
    sys.saveState(s);
    std::cout<<s.str()<<std::endl;
    sys.loadState(s);
    
    BOOST_CHECK( boost::num_vertices(sys.m_cluster) == 5 );
    BOOST_CHECK( boost::num_edges(sys.m_cluster) == 0 );
    BOOST_CHECK( sys.m_cluster.numClusters() == 0 );
    
    std::cout<<"load state done"<<std::endl;    
}

BOOST_AUTO_TEST_SUITE_END();
