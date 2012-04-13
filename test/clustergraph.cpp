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

    Graph g2 = g1.createCluster();
    BOOST_CHECK(!g2.isRoot());
    BOOST_CHECK(g1==g2.parent());

    Graph g3 = g2.createCluster();
    BOOST_CHECK(g1==g3.root());
    BOOST_CHECK(g3.numClusters() == 0);
    BOOST_CHECK(g3.clusters().first == g3.clusters().second);

    Graph g4 = g2.createCluster();
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

BOOST_AUTO_TEST_CASE(object_property_handling) {

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
    
}

BOOST_AUTO_TEST_CASE(movevertex) {

    //check vertex reclustering
    Graph g;
    fusion::vector<LocalVertex, GlobalVertex> v1 = g.addVertex();
    fusion::vector<LocalVertex, GlobalVertex> v2 = g.addVertex();
    fusion::vector<LocalVertex, GlobalVertex> v3 = g.addVertex();

    //g.addEdge(v1, v2);
    //g.addEdge(v2,v3);

}
