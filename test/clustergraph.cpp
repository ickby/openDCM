/*
    <one line to give the program's name and a brief idea of what it does.>
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

#include "clustergraph.hpp"
#include "property.hpp"

#include <boost/fusion/sequence.hpp>
#include <boost/mpl/vector.hpp>
#include <utility>

  #include <iostream>

//  #include <boost/fusion/container.hpp>
//  #include <boost/fusion/sequence.hpp>
//  #include <boost/smart_ptr/shared_ptr.hpp>

#define BOOST_TEST_MODULE ClusterGraph
#include <boost/test/unit_test.hpp>

using namespace dcm;
namespace mpl = boost::mpl;

struct test_edge_property {
    typedef dcm::edge_property_tag kind;
    typedef int type;
};

struct test_vertex_property {
    typedef dcm::vertex_property_tag kind;
    typedef int type;
};

struct test_object1 {
    int value;
};
struct test_object2 {
    int value;
};

typedef ClusterGraph<mpl::vector<test_edge_property>,
mpl::vector<test_vertex_property>, mpl::vector<test_object1, test_object2> > Graph;

BOOST_AUTO_TEST_CASE(subclustering) {

    Graph g1;
    BOOST_CHECK( g1.isRoot() );

    Graph g2 = g1.createCluster().first;
    BOOST_CHECK(!g2.isRoot());
    BOOST_CHECK(g1==g2.parent());

    Graph g3 = g2.createCluster().first;
    BOOST_CHECK(g1==g3.root());
    BOOST_CHECK(g3.numClusters() == 0);
    BOOST_CHECK(g3.clusters().first == g3.clusters().second);

    Graph g4 = g2.createCluster().first;
    Graph::cluster_iterator it,end;
    boost::tie(it,end) = g2.clusters();
    BOOST_CHECK(it != end);
    BOOST_CHECK(it->second->parent() == (it++)->second->parent());
    BOOST_CHECK(g2.numClusters() == 2);
}

BOOST_AUTO_TEST_CASE(creation_handling) {
  
  Graph g1;
  fusion::vector<LocalVertex, GlobalVertex> res1 = g1.addVertex();
  fusion::vector<LocalVertex, GlobalVertex> res2 = g1.addVertex();
  BOOST_CHECK(fusion::at_c<0>(res1) != fusion::at_c<0>(res2));
  BOOST_CHECK(fusion::at_c<1>(res1) != fusion::at_c<1>(res2));
  std::pair<LocalVertex, bool> loc = g1.getLocalVertex( fusion::at_c<1>(res1) );
  BOOST_CHECK(loc.second);
  BOOST_CHECK(fusion::at_c<0>(res1) == loc.first);
  BOOST_CHECK(fusion::at_c<1>(res1) == g1.getGlobalVertex(fusion::at_c<0>(res1)));
  
  fusion::vector<LocalEdge, GlobalEdge, bool> edge1 = g1.addEdge( fusion::at_c<0>(res1), fusion::at_c<0>(res2) );
  fusion::vector<LocalEdge, GlobalEdge, bool> edge2 = g1.addEdge( fusion::at_c<0>(res1), fusion::at_c<0>(res2) );
  BOOST_CHECK(fusion::at_c<2>(edge1));
  std::pair<LocalEdge, bool> loc2 = g1.getLocalEdge(fusion::at_c<1>(edge1));
  BOOST_CHECK(loc2.second);
  BOOST_CHECK(fusion::at_c<0>(edge1) == loc2.first);
  BOOST_CHECK(fusion::at_c<1>(edge1) == *(g1.getGlobalEdges(fusion::at_c<0>(edge1)).first));
  
  BOOST_CHECK(fusion::at_c<2>(edge2));
  BOOST_CHECK(fusion::at_c<0>(edge2) == fusion::at_c<0>(edge1));
  BOOST_CHECK(fusion::at_c<1>(edge2) == fusion::at_c<1>(edge1));
};

BOOST_AUTO_TEST_CASE(object_handling) {
    
    Graph g1;
    fusion::vector<LocalVertex, GlobalVertex> v1c = g1.addVertex();
    fusion::vector<LocalVertex, GlobalVertex> v2c = g1.addVertex();
    LocalVertex v1 = fusion::at_c<0>(v1c);
    LocalVertex v2 = fusion::at_c<0>(v2c);
    
    fusion::vector<LocalEdge, GlobalEdge, bool> seq = g1.addEdge(v1,v2);
    LocalEdge e1 = fusion::at_c<0>(seq);
    bool inserted = fusion::at_c<2>(seq);
    BOOST_CHECK(inserted);
    
    boost::shared_ptr<test_object1> o1( new test_object1 );
    boost::shared_ptr<test_object2> o2( new test_object2 );
    o1->value = 1;

    BOOST_CHECK( !g1.getObject<test_object1>(v1) );
    BOOST_CHECK( !g1.getObject<test_object2>(v1) );
    
    g1.setObject(v1, o1);    
    BOOST_CHECK( g1.getObject<test_object1>(v1) );
    BOOST_CHECK( !g1.getObject<test_object2>(v1));    
    BOOST_CHECK( o1 == g1.getObject<test_object1>(v1));
    BOOST_CHECK( o1->value == g1.getObject<test_object1>(v1)->value);
    
    g1.setObject(v1,o2);        
    BOOST_CHECK( !g1.getObject<test_object1>(v1) );
    BOOST_CHECK( g1.getObject<test_object2>(v1));
    BOOST_CHECK( o2 == g1.getObject<test_object2>(v1));
    

    BOOST_CHECK( !g1.getObject<test_object1>(e1) );
    BOOST_CHECK( !g1.getObject<test_object2>(e1) );
    
    g1.setObject(e1, o1);
    BOOST_CHECK( g1.getObject<test_object1>(e1) );
    BOOST_CHECK( !g1.getObject<test_object2>(e1));    
    BOOST_CHECK( o1 == g1.getObject<test_object1>(e1));
    BOOST_CHECK( o1->value == g1.getObject<test_object1>(e1)->value);

    g1.setObject(e1,o2);
    BOOST_CHECK( !g1.getObject<test_object1>(e1) );
    BOOST_CHECK( g1.getObject<test_object2>(e1));
    BOOST_CHECK( o2 == g1.getObject<test_object2>(e1));
    
    BOOST_CHECK( o2 == *g1.getObjects<test_object2>(e1).first );
    BOOST_CHECK( *g1.getObjects<test_object2>(e1).first != *g1.getObjects<test_object2>(e1).second );
    BOOST_CHECK( *(++g1.getObjects<test_object2>(e1).first) == *g1.getObjects<test_object2>(e1).second );
}

BOOST_AUTO_TEST_CASE(property_handling) {
    
    Graph g1;
    fusion::vector<LocalVertex, GlobalVertex> v1c = g1.addVertex();
    fusion::vector<LocalVertex, GlobalVertex> v2c = g1.addVertex();
    LocalVertex v1 = fusion::at_c<0>(v1c);
    LocalVertex v2 = fusion::at_c<0>(v2c);
    
    fusion::vector<LocalEdge, GlobalEdge, bool> seq = g1.addEdge(v1,v2);
    GlobalEdge e1 = fusion::at_c<1>(seq);
    bool inserted = fusion::at_c<2>(seq);
    BOOST_CHECK(inserted);
    

    g1.setProperty<test_vertex_property>(v1, 7);    
    BOOST_CHECK( g1.getProperty<test_vertex_property>(v1) == 7 );
 
    g1.setProperty<test_edge_property>(e1, 3);    
    BOOST_CHECK( g1.getProperty<test_edge_property>(e1) == 3 );
    BOOST_CHECK( g1.getProperty<test_edge_property>(g1.getLocalEdge(e1).first) == 3 );    
    
    LocalEdge e = fusion::at_c<0>(seq);
    BOOST_CHECK( 3 == *g1.getProperties<test_edge_property>(e).first );
    BOOST_CHECK( *g1.getProperties<test_vertex_property>(e).first != *g1.getProperties<test_vertex_property>(e).second );
    BOOST_CHECK( *(++g1.getProperties<test_vertex_property>(e).first) == *g1.getProperties<test_vertex_property>(e).second );
}

BOOST_AUTO_TEST_CASE(vertex_to_subcluster) {

    //check vertex reclustering
    Graph g;
    fusion::vector<LocalVertex, GlobalVertex> v1 = g.addVertex();
    fusion::vector<LocalVertex, GlobalVertex> v2 = g.addVertex();
    fusion::vector<LocalVertex, GlobalVertex> v3 = g.addVertex();
    
    std::pair<Graph&, LocalVertex> res = g.createCluster();
    fusion::vector<LocalVertex, GlobalVertex> v4 = res.first.addVertex();

    GlobalEdge e1 = fusion::at_c<1>( g.addEdge(fusion::at_c<0>(v1),fusion::at_c<0>(v2)) );
    GlobalEdge e2 = fusion::at_c<1>( g.addEdge(fusion::at_c<0>(v2),fusion::at_c<0>(v3)) );
    GlobalEdge e3 = fusion::at_c<1>( g.addEdge(fusion::at_c<0>(v3),fusion::at_c<0>(v1)) );
    
    GlobalEdge e4 = fusion::at_c<1>( g.addEdge(fusion::at_c<1>(v1),fusion::at_c<1>(v4)) );
    GlobalEdge e5 = fusion::at_c<1>( g.addEdge(fusion::at_c<1>(v3),fusion::at_c<1>(v4)) );
    
    g.setProperty<test_vertex_property>(fusion::at_c<0>(v1), 25);
    g.setProperty<test_edge_property>(e1, 1);
    g.setProperty<test_edge_property>(e2, 2);
    g.setProperty<test_edge_property>(e3, 3);
    g.setProperty<test_edge_property>(e4, 4);
    g.setProperty<test_edge_property>(e5, 5);
    
    //there should be edges between cluster and v1 v3
    BOOST_CHECK( g.edge(fusion::at_c<0>(v1), res.second).second );
    BOOST_CHECK( g.edge(fusion::at_c<0>(v3), res.second).second );
    BOOST_CHECK( *g.getProperties<test_edge_property>( g.edge(fusion::at_c<0>(v1), res.second).first).first == 4 );
    BOOST_CHECK( *g.getProperties<test_edge_property>( g.edge(fusion::at_c<0>(v3), res.second).first).first == 5 ); 
    
    //inside subcluster vertex v4 should not have any edges
    BOOST_CHECK( boost::out_degree(fusion::at_c<0>(v4), res.first) == 0 );
    
    LocalVertex nv = g.vertexToSubcluster(fusion::at_c<0>(v1), res.second);
    
    BOOST_CHECK( res.first.getProperty<test_vertex_property>(nv) == 25 );
    BOOST_CHECK( res.first.getGlobalVertex(nv) == fusion::at_c<1>(v1) );
    BOOST_CHECK( boost::out_degree(fusion::at_c<0>(v4), res.first) == 1 );
    
    //in the local edge between new vertex and v4 should be one global edge
    LocalEdge et0 = res.first.edge(fusion::at_c<0>(v4), nv).first;
    typedef typename Graph::property_iterator<test_edge_property> prop_iter;
    std::pair< prop_iter, prop_iter > it = res.first.getProperties<test_edge_property>(et0);
    BOOST_CHECK( *it.first == 4);
    BOOST_CHECK( ++it.first == it.second );
    
    //v3 should have two global edges to the cluster, his own and the one that was to v1 before
    LocalEdge et1 = g.edge(fusion::at_c<0>(v3), res.second).first;
    it = g.getProperties<test_edge_property>(et1);
    BOOST_CHECK( *it.first == 5);
    BOOST_CHECK( *(++it.first) == 3);
    BOOST_CHECK( ++it.first == it.second );
    
    //v2 should have one global edge to the cluster, the one that was to v1 before
    LocalEdge et2 = g.edge(fusion::at_c<0>(v2), res.second).first;
    it = g.getProperties<test_edge_property>(et2);
    BOOST_CHECK( *it.first == 1);
    BOOST_CHECK( ++it.first == it.second );

}
