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

#ifdef _WIN32
	//warning about to long decoraded names, won't affect the code correctness
	#pragma warning( disable : 4503 )
#endif

#include "opendcm/core/clustergraph.hpp"
#include "opendcm/core/filtergraph.hpp"

#include <boost/graph/undirected_dfs.hpp>

#include <boost/test/unit_test.hpp>

using namespace dcm;
using namespace graph;
namespace mpl = boost::mpl;

BOOST_AUTO_TEST_SUITE(ClusterGraph_test_suit);

struct test_edge_property {
    typedef int type;
    struct default_value {
	int operator()() {
	  return 2;
	};
    };
    struct change_tracking {};
};

struct test_edge_property2 {
    typedef int type;
};

struct test_globaledge_property {
    typedef int type;
    struct change_tracking {};
};

struct test_vertex_property {
    typedef int type;
    struct change_tracking {};
};

struct test_cluster_property {
    typedef int type;
};

typedef ClusterGraph<mpl::vector2<test_edge_property, test_edge_property2>, mpl::vector1<test_globaledge_property>,
        mpl::vector1<test_vertex_property>, mpl::vector1<test_cluster_property> > Graph;

struct delete_functor {
    std::stringstream stream;
    
    void operator()(GlobalVertex v) {
      stream<<"v";
    }
    void operator()(GlobalEdge e) {
      stream<<"e";
    }
    void operator()(std::shared_ptr<Graph> v) {
      stream<<"c";
    }
};

BOOST_AUTO_TEST_CASE(subclustering) {
	
    std::shared_ptr<Graph> g1 = std::shared_ptr<Graph>(new Graph);
    BOOST_CHECK(g1->isRoot());

    std::pair<  std::shared_ptr<Graph>, LocalVertex> sub = g1->createCluster();
    std::shared_ptr<Graph> g2 = sub.first;
    BOOST_CHECK(!g2->isRoot());
    BOOST_CHECK(g1==g2->parent());
    BOOST_CHECK(g1->numClusters() == 1);
    BOOST_CHECK(g2 == g1->getVertexCluster(sub.second));
    BOOST_CHECK(sub.second == g1->getClusterVertex(sub.first));
	
    std::shared_ptr<Graph> g3 = g2->createCluster().first;
    BOOST_CHECK(g1==g3->root());
    BOOST_CHECK(g3->numClusters() == 0);
    BOOST_CHECK(g3->clusters().first == g3->clusters().second);
    
    g2->addVertex();

    std::shared_ptr<Graph> g4 = g2->createCluster().first;
    Graph::cluster_iterator it,end;
    boost::tie(it,end) = g2->clusters();
    BOOST_CHECK(it != end);
    BOOST_CHECK(it->second->parent() == (it++)->second->parent());
    BOOST_CHECK(g2->numClusters() == 2);
    BOOST_CHECK(g2->vertexCount()==3);
    
    delete_functor f;
    g2->removeCluster(g3, f);
    boost::tie(it,end) = g2->clusters();
    BOOST_CHECK(it != end);
    BOOST_CHECK(g2->numClusters() == 1);
    BOOST_CHECK(g2->vertexCount()==2);
    BOOST_CHECK( !f.stream.str().compare("c") );
    
    delete_functor f2;
    g1->removeCluster(sub.first, f2);
    BOOST_CHECK( !f2.stream.str().compare("cvc") );
    

}

BOOST_AUTO_TEST_CASE(creation_handling) {

    std::shared_ptr<Graph> g1 = std::shared_ptr<Graph>(new Graph);
    fusion::vector<LocalVertex, GlobalVertex> sub1 = g1->addVertex();
    fusion::vector<LocalVertex, GlobalVertex> sub2 = g1->addVertex();
    BOOST_CHECK(fusion::at_c<0>(sub1) != fusion::at_c<0>(sub2));
    BOOST_CHECK(fusion::at_c<1>(sub1) != fusion::at_c<1>(sub2));
    std::pair<LocalVertex, bool> loc = g1->getLocalVertex(fusion::at_c<1>(sub1));
    BOOST_CHECK(loc.second);
    BOOST_CHECK(fusion::at_c<0>(sub1) == loc.first);
    BOOST_CHECK(fusion::at_c<1>(sub1) == g1->getGlobalVertex(fusion::at_c<0>(sub1)));

    fusion::vector<LocalEdge, GlobalEdge, bool> edge1 = g1->addEdge(fusion::at_c<0>(sub1), fusion::at_c<0>(sub2));
    fusion::vector<LocalEdge, GlobalEdge, bool> edge2 = g1->addEdge(fusion::at_c<0>(sub1), fusion::at_c<0>(sub2));
    BOOST_CHECK(fusion::at_c<2>(edge1));
    std::pair<LocalEdge, bool> loc2 = g1->getLocalEdge(fusion::at_c<1>(edge1));
    BOOST_CHECK(loc2.second);
    BOOST_CHECK(fusion::at_c<0>(edge1) == loc2.first);
    BOOST_CHECK(fusion::at_c<1>(edge1) == *(g1->getGlobalEdges(fusion::at_c<0>(edge1)).first));

    BOOST_CHECK(fusion::at_c<2>(edge2));
    BOOST_CHECK(fusion::at_c<0>(edge2) == fusion::at_c<0>(edge1));
    BOOST_CHECK(fusion::at_c<1>(edge2) != fusion::at_c<1>(edge1));
    
    //check edge creation when 1 vertex is in a cluster
    std::pair<std::shared_ptr<Graph>, LocalVertex> nc = g1->createCluster();
    fusion::vector<LocalVertex, GlobalVertex> sub3 = nc.first->addVertex();
    fusion::vector<LocalEdge, GlobalEdge, bool> edge3 = g1->addEdge(fusion::at_c<1>(sub2), fusion::at_c<1>(sub3));
    BOOST_CHECK( fusion::at_c<2>(edge3) );
    BOOST_CHECK( boost::source(fusion::at_c<0>(edge3), g1->getDirectAccess()) == fusion::at_c<0>(sub2) );
    BOOST_CHECK( boost::target(fusion::at_c<0>(edge3), g1->getDirectAccess()) == nc.second );
    
};

BOOST_AUTO_TEST_CASE(removing) {

    std::shared_ptr<Graph> g1 = std::shared_ptr<Graph>(new Graph);
    fusion::vector<LocalVertex, GlobalVertex> res1 = g1->addVertex();
    fusion::vector<LocalVertex, GlobalVertex> res2 = g1->addVertex();
    fusion::vector<LocalVertex, GlobalVertex> res3 = g1->addVertex();
    
    fusion::vector<LocalEdge, GlobalEdge, bool> res4 = g1->addEdge(fusion::at_c<0>(res1), fusion::at_c<0>(res2));
    fusion::vector<LocalEdge, GlobalEdge, bool> res5 = g1->addEdge(fusion::at_c<0>(res2), fusion::at_c<0>(res3));
    
    g1->removeVertex(fusion::at_c<0>(res1));    
    BOOST_CHECK( g1->edgeCount() == 1 );
    BOOST_CHECK( g1->vertexCount() == 2 );
    
    g1->removeVertex(fusion::at_c<1>(res2));    
    BOOST_CHECK( g1->edgeCount() == 0 );
    BOOST_CHECK( g1->vertexCount() == 1 );
    
    //create subcluster with two vertices
    std::shared_ptr<Graph> g2 = g1->createCluster().first;
    res1 = g2->addVertex();
    res2 = g2->addVertex();
    
    //connect toplevel vertex with both subcluster vertices, so that there should be one local edge in top 
    //cluster with to global ones inside
    res4 = g1->addEdge(fusion::at_c<1>(res3), fusion::at_c<1>(res1));
    res5 = g1->addEdge(fusion::at_c<1>(res3), fusion::at_c<1>(res2));
    BOOST_CHECK( g1->edgeCount() == 1 );
    
    g1->removeVertex(fusion::at_c<1>(res1));
    BOOST_CHECK( g1->vertexCount() == 2 ); //one normal one cluster
    BOOST_CHECK( g2->vertexCount() == 1 ); //one normal one cluster
    BOOST_CHECK( g1->edgeCount() == 1 );
    
    g1->removeVertex(fusion::at_c<1>(res2));
    BOOST_CHECK( g1->vertexCount() == 2 ); //one normal one cluster
    BOOST_CHECK( g2->vertexCount() == 0 ); //one normal one cluster
    BOOST_CHECK( g1->edgeCount() == 0 );
    
    res1 = g2->addVertex();
    res2 = g2->addVertex();
    res4 = g1->addEdge(fusion::at_c<1>(res3), fusion::at_c<1>(res1));
    res5 = g1->addEdge(fusion::at_c<1>(res3), fusion::at_c<1>(res2));
    fusion::vector<LocalEdge, GlobalEdge, bool> res6 = g2->addEdge(fusion::at_c<0>(res1), fusion::at_c<0>(res2));
    
    BOOST_CHECK( g2->edgeCount() == 1 );
    BOOST_CHECK( g1->edgeCount() == 1 );
    BOOST_CHECK( g1->getGlobalEdgeCount(fusion::at_c<0>(res4)) == 2);
    
    g1->removeEdge(fusion::at_c<1>(res4));
    BOOST_CHECK( g2->edgeCount() == 1 );
    BOOST_CHECK( g1->edgeCount() == 1 );
    BOOST_CHECK( g1->getGlobalEdgeCount(fusion::at_c<0>(res4)) == 1);
    
    g1->removeEdge(fusion::at_c<1>(res5));
    BOOST_CHECK( g2->edgeCount() == 1 );
    BOOST_CHECK( g1->edgeCount() == 0 );
    
};

template<typename T>
void pretty(T t) {
    std::cout<<__PRETTY_FUNCTION__<<std::endl;
};

template<typename P, typename G, typename K>
void prettymap(property_map<P,G,K> map) {
    std::cout<<__PRETTY_FUNCTION__<<std::endl;
}    

BOOST_AUTO_TEST_CASE(property_handling) {

    std::shared_ptr<Graph> g1 = std::shared_ptr<Graph>(new Graph);
    g1->setProperty<test_cluster_property>(5);
    BOOST_CHECK(g1->getProperty<test_cluster_property>() == 5);


    fusion::vector<LocalVertex, GlobalVertex> v1c = g1->addVertex();
    fusion::vector<LocalVertex, GlobalVertex> v2c = g1->addVertex();
    LocalVertex v1 = fusion::at_c<0>(v1c);
    LocalVertex v2 = fusion::at_c<0>(v2c);

    fusion::vector<LocalEdge, GlobalEdge, bool> seq = g1->addEdge(v1,v2);
    LocalEdge e = fusion::at_c<0>(seq);
    GlobalEdge e1 = fusion::at_c<1>(seq);
    bool inserted = fusion::at_c<2>(seq);
    BOOST_CHECK(inserted);


    BOOST_CHECK(!g1->isPropertyChanged<test_vertex_property>(v1));
    g1->setProperty<test_vertex_property>(v1, 7);
    BOOST_CHECK(g1->isPropertyChanged<test_vertex_property>(v1));
    BOOST_CHECK(g1->getProperty<test_vertex_property>(v1) == 7);
    g1->acknowledgePropertyChange<test_vertex_property>(v1);
    BOOST_CHECK(!g1->isPropertyChanged<test_vertex_property>(v1));

    g1->setProperty<test_vertex_property>(fusion::at_c<1>(v1c), 2);
    BOOST_CHECK(g1->getProperty<test_vertex_property>(v1) == 2);

    BOOST_CHECK(g1->getProperty<test_edge_property>(e) == 2);
    g1->setProperty<test_edge_property>(e, 3);
    BOOST_CHECK(g1->getProperty<test_edge_property>(e) == 3);
    
    g1->setProperty<test_globaledge_property>(e1, 5);
    BOOST_CHECK(g1->getProperty<test_globaledge_property>(e1) == 5);

    g1->setProperty<test_edge_property2>(e, 2);
    BOOST_CHECK(g1->getProperty<test_edge_property>(e) == 3);
    BOOST_CHECK(g1->getProperty<test_edge_property2>(e) == 2);

    //test property maps
    property_map<test_edge_property2, Graph> pmap(g1);
    BOOST_CONCEPT_ASSERT(( boost::WritablePropertyMapConcept<property_map<test_edge_property2, Graph>, graph::LocalEdge> ));
    BOOST_CONCEPT_ASSERT(( boost::ReadWritePropertyMapConcept<property_map<test_edge_property2, Graph>, graph::LocalEdge> ));

    BOOST_CHECK(boost::get(pmap, e) == 2);
    boost::put(pmap, e, 7);
    BOOST_CHECK(boost::get(pmap, e) == 7);
    BOOST_CHECK(boost::at(pmap, e) == 7);
    boost::at(pmap, e) = 5;
    BOOST_CHECK(boost::get(pmap, e) == 5);
    
    //see it the default properties work in an boost algorithm
    property_map<Index, Graph, LocalVertex> imap(g1);
    property_map<Color, Graph, LocalEdge>   ie_cpm(g1);
    property_map<Color, Graph, LocalVertex> iv_cpm(g1);
    
    boost::undirected_dfs(g1->getDirectAccess(), boost::dfs_visitor<boost::null_visitor>(), iv_cpm, ie_cpm, *g1->vertices().first);

}

BOOST_AUTO_TEST_CASE(property_iterating) {

    std::shared_ptr<Graph> g1 = std::shared_ptr<Graph>(new Graph);
 
    fusion::vector<LocalVertex, GlobalVertex> v1c = g1->addVertex();
    fusion::vector<LocalVertex, GlobalVertex> v2c = g1->addVertex();
    fusion::vector<LocalVertex, GlobalVertex> v3c = g1->addVertex();

    fusion::vector<LocalEdge, GlobalEdge, bool> e1c = g1->addEdge(fusion::at_c<0>(v1c),fusion::at_c<0>(v2c));
    fusion::vector<LocalEdge, GlobalEdge, bool> e2c = g1->addEdge(fusion::at_c<0>(v2c),fusion::at_c<0>(v3c));
    fusion::vector<LocalEdge, GlobalEdge, bool> e3c = g1->addEdge(fusion::at_c<0>(v3c),fusion::at_c<0>(v1c));

    //at the begining every edge and vertex should be marked as changed as they are newly created
    typedef boost::filter_iterator<Graph::property_changes, Graph::local_edge_iterator> eit1;
    std::pair<eit1, eit1> p = g1->filterRange<Graph::property_changes>(boost::edges(g1->getDirectAccess()));
    int c = 0;
    for(; p.first != p.second; ++p.first) {
        BOOST_CHECK(g1->edgeChanged(*p.first));
        g1->acknowledgeEdgeChanges(*p.first);
        ++c;
    }
    
    BOOST_CHECK( c == 3 );
    
    typedef boost::filter_iterator<Graph::property_changes, Graph::local_vertex_iterator> vit1;
    std::pair<vit1, vit1> v = g1->filterRange<Graph::property_changes>(boost::vertices(g1->getDirectAccess()));
    c = 0;
    for(; v.first != v.second; ++v.first) {
        g1->acknowledgePropertyChanges(*v.first);
        ++c;
    }
    
    BOOST_CHECK( c == 3 );
    
    //We acknowledged all changes, lets see if this worked
    p = g1->filterRange<Graph::property_changes>(boost::edges(g1->getDirectAccess()));
    BOOST_CHECK(p.first == p.second);    
    v = g1->filterRange<Graph::property_changes>(boost::vertices(g1->getDirectAccess()));
    BOOST_CHECK(v.first == v.second);
    
    //see if we track changes correctly
    g1->setProperty<test_vertex_property>(fusion::at_c<1>(v2c), 5);
    
    v = g1->filterRange<Graph::property_changes>(boost::vertices(g1->getDirectAccess()));
    BOOST_CHECK( *v.first == fusion::at_c<0>(v2c) );
    BOOST_CHECK( ++v.first == v.second );
    
    typedef boost::filter_iterator<Graph::property_changed<test_vertex_property>, Graph::local_vertex_iterator> vit2;
    std::pair<vit2, vit2> v2 = g1->filterRange<Graph::property_changed<test_vertex_property> >(boost::vertices(g1->getDirectAccess()));
    BOOST_CHECK( *v2.first == fusion::at_c<0>(v2c) );
    BOOST_CHECK( ++v2.first == v2.second );
    
    typedef boost::filter_iterator<Graph::property_value<test_vertex_property, 5>, Graph::local_vertex_iterator> vit3;
    std::pair<vit3, vit3> v3 = g1->filterRange<Graph::property_value<test_vertex_property, 5> >(boost::vertices(g1->getDirectAccess()));
    BOOST_CHECK( *v3.first == fusion::at_c<0>(v2c) );
    BOOST_CHECK( ++v3.first == v3.second );

}

BOOST_AUTO_TEST_CASE(move_vertex) {

    //check vertex reclustering
    std::shared_ptr<Graph> g = std::shared_ptr<Graph>(new Graph);
    fusion::vector<LocalVertex, GlobalVertex> v1 = g->addVertex();
    fusion::vector<LocalVertex, GlobalVertex> v2 = g->addVertex();
    fusion::vector<LocalVertex, GlobalVertex> v3 = g->addVertex();

    std::pair<std::shared_ptr<Graph>, LocalVertex> sub = g->createCluster();
    fusion::vector<LocalVertex, GlobalVertex> v4 = sub.first->addVertex();

    GlobalEdge e1 = fusion::at_c<1>(g->addEdge(fusion::at_c<0>(v1),fusion::at_c<0>(v2)));
    GlobalEdge e2 = fusion::at_c<1>(g->addEdge(fusion::at_c<0>(v2),fusion::at_c<0>(v3)));
    GlobalEdge e3 = fusion::at_c<1>(g->addEdge(fusion::at_c<0>(v3),fusion::at_c<0>(v1)));

    GlobalEdge e4 = fusion::at_c<1>(g->addEdge(fusion::at_c<1>(v1),fusion::at_c<1>(v4)));
    GlobalEdge e5 = fusion::at_c<1>(g->addEdge(fusion::at_c<1>(v3),fusion::at_c<1>(v4)));

    g->setProperty<test_vertex_property>(fusion::at_c<0>(v1),25);

    //there should be edges between cluster and v1 v3
    BOOST_CHECK(g->edge(fusion::at_c<0>(v1), sub.second).second);
    BOOST_CHECK(g->edge(fusion::at_c<0>(v3), sub.second).second);
    BOOST_CHECK(*g->getGlobalEdges(g->edge(fusion::at_c<0>(v1), sub.second).first).first == e4);
    BOOST_CHECK(*g->getGlobalEdges(g->edge(fusion::at_c<0>(v3), sub.second).first).first == e5);

    //inside subcluster vertex v4 should not have any edges
    BOOST_CHECK(sub.first->vertexCount() == 1);
    BOOST_CHECK(sub.first->edgeCount() == 0);
    BOOST_CHECK(sub.first->outDegree(fusion::at_c<0>(v4)) == 0);

    LocalVertex nv = g->moveToSubcluster(fusion::at_c<0>(v1), sub.second);

    BOOST_CHECK(sub.first->getProperty<test_vertex_property>(nv) == 25);
    BOOST_CHECK(sub.first->getGlobalVertex(nv) == fusion::at_c<1>(v1));
    BOOST_CHECK(sub.first->vertexCount() == 2);
    BOOST_CHECK(sub.first->edgeCount() == 1);
    BOOST_CHECK(sub.first->outDegree(fusion::at_c<0>(v4)) == 1);
    BOOST_CHECK(g->vertexCount() == 3);
    BOOST_CHECK(g->edgeCount() == 3);

    //in the local edge between new vertex and v4 should be one global edge
    LocalEdge et0;
    bool is_edge;
    boost::tie(et0, is_edge) = sub.first->edge(fusion::at_c<0>(v4), nv);
    BOOST_CHECK(is_edge);
    std::pair< Graph::global_edge_iterator, Graph::global_edge_iterator > it = sub.first->getGlobalEdges(et0);
    BOOST_CHECK(*it.first == e4);
    BOOST_CHECK(++it.first == it.second);

    //v3 should have two global edges to the cluster, his own and the one that was to v1 before
    LocalEdge et1 = g->edge(fusion::at_c<0>(v3), sub.second).first;
    it = g->getGlobalEdges(et1);
    BOOST_CHECK(*it.first == e5);
    BOOST_CHECK(*(++it.first) == e3);
    BOOST_CHECK(++it.first == it.second);

    //v2 should have one global edge to the cluster, the one that was to v1 before
    LocalEdge et2 = g->edge(fusion::at_c<0>(v2), sub.second).first;
    it = g->getGlobalEdges(et2);
    BOOST_CHECK(*it.first == e1);
    BOOST_CHECK(++it.first == it.second);


    LocalVertex nv1 = sub.first->moveToParent(nv);

    //everything should be like in the beginning, first check subcluster
    BOOST_CHECK(sub.first->outDegree(fusion::at_c<0>(v4)) == 0);
    BOOST_CHECK(sub.first->vertexCount() == 1);

    BOOST_CHECK(g->vertexCount() == 4);
    BOOST_CHECK(g->edgeCount() == 5);
    BOOST_CHECK(g->outDegree(sub.second) == 2);
    BOOST_CHECK(g->outDegree(nv1) == 3);
    BOOST_CHECK(g->edge(nv1, sub.second).second);
    BOOST_CHECK(!g->edge(fusion::at_c<0>(v2), sub.second).second);
    BOOST_CHECK(g->edge(fusion::at_c<0>(v3), sub.second).second);

    //nv1 to cluster should have one global edge
    LocalEdge et3 = g->edge(nv1, sub.second).first;
    it = g->getGlobalEdges(et3);
    BOOST_CHECK(*it.first == e4);
    BOOST_CHECK(++it.first == it.second);

    //nv1 to v2 should have one global edge
    LocalEdge et4 = g->edge(nv1, fusion::at_c<0>(v2)).first;
    it = g->getGlobalEdges(et4);
    BOOST_CHECK(*it.first == e1);
    BOOST_CHECK(++it.first == it.second);

    //nv1 to v3 should have one global edge
    LocalEdge et5 = g->edge(nv1, fusion::at_c<0>(v3)).first;
    it = g->getGlobalEdges(et5);
    BOOST_CHECK(*it.first == e3);
    BOOST_CHECK(++it.first == it.second);

    //nv2 to v3 should have one global edge
    LocalEdge et6 = g->edge(fusion::at_c<0>(v2), fusion::at_c<0>(v3)).first;
    it = g->getGlobalEdges(et6);
    BOOST_CHECK(*it.first == e2);
    BOOST_CHECK(++it.first == it.second);

    //nv3 to cluster should have one global edge
    LocalEdge et7 = g->edge(fusion::at_c<0>(v3), sub.second).first;
    it = g->getGlobalEdges(et7);
    BOOST_CHECK(*it.first == e5);
    BOOST_CHECK(++it.first == it.second);
}

BOOST_AUTO_TEST_CASE(move_cluster) {

    //check vertex reclustering
    std::shared_ptr<Graph> g = std::shared_ptr<Graph>(new Graph);
    fusion::vector<LocalVertex, GlobalVertex> v1 = g->addVertex();
    fusion::vector<LocalVertex, GlobalVertex> v2 = g->addVertex();

    std::pair<std::shared_ptr<Graph>, LocalVertex> sub1 = g->createCluster();
    fusion::vector<LocalVertex, GlobalVertex> v3 = sub1.first->addVertex();
    std::pair<std::shared_ptr<Graph>, LocalVertex> sub2 = g->createCluster();
    fusion::vector<LocalVertex, GlobalVertex> v4 = sub2.first->addVertex();
    std::pair<std::shared_ptr<Graph>, LocalVertex> sub3 = g->createCluster();
    fusion::vector<LocalVertex, GlobalVertex> v5 = sub3.first->addVertex();

    GlobalEdge e12 = fusion::at_c<1>(g->addEdge(fusion::at_c<0>(v1),fusion::at_c<0>(v2)));
    GlobalEdge e13 = fusion::at_c<1>(g->addEdge(fusion::at_c<1>(v1),fusion::at_c<1>(v3)));
    GlobalEdge e23 = fusion::at_c<1>(g->addEdge(fusion::at_c<1>(v2),fusion::at_c<1>(v3)));
    GlobalEdge e24 = fusion::at_c<1>(g->addEdge(fusion::at_c<1>(v2),fusion::at_c<1>(v4)));
    GlobalEdge e34 = fusion::at_c<1>(g->addEdge(fusion::at_c<1>(v4),fusion::at_c<1>(v3)));
    GlobalEdge e45 = fusion::at_c<1>(g->addEdge(fusion::at_c<1>(v4),fusion::at_c<1>(v5)));
    GlobalEdge e35 = fusion::at_c<1>(g->addEdge(fusion::at_c<1>(v3),fusion::at_c<1>(v5)));

    LocalVertex nv = g->moveToSubcluster(sub2.second, sub3.second);

    //first check local cluster
    BOOST_CHECK(g->edgeCount() == 5);
    BOOST_CHECK(g->vertexCount() == 4);

    //sub1 to sub2 should have two global edge's
    LocalEdge et1 = g->edge(sub1.second, sub3.second).first;
    std::pair< Graph::global_edge_iterator, Graph::global_edge_iterator > it;
    it = g->getGlobalEdges(et1);
    BOOST_CHECK(*it.first == e35);
    BOOST_CHECK(*(++it.first) == e34);
    BOOST_CHECK(++it.first == it.second);

    //v2 to sub3 should have one global edge
    LocalEdge et2 = g->edge(fusion::at_c<0>(v2), sub3.second).first;
    it = g->getGlobalEdges(et2);
    BOOST_CHECK(*it.first == e24);
    BOOST_CHECK(++it.first == it.second);

    //check subcluster 3
    BOOST_CHECK(sub3.first->edgeCount() == 1);
    BOOST_CHECK(sub3.first->vertexCount() == 2);

    //nv to v5 should have one global edge
    LocalEdge et3 = sub3.first->edge(nv, fusion::at_c<0>(v5)).first;
    it = sub3.first->getGlobalEdges(et3);
    BOOST_CHECK(*it.first == e45);
    BOOST_CHECK(++it.first == it.second);


    LocalVertex nc = sub3.first->moveToParent(nv);

    //everything need to be like in initial state
    BOOST_CHECK(g->edgeCount() == 7);
    BOOST_CHECK(g->vertexCount() == 5);
    BOOST_CHECK(sub3.first->edgeCount() == 0);
    BOOST_CHECK(sub3.first->vertexCount() == 1);

    //nc to sub1 should have one global edge
    std::pair<LocalEdge, bool> res = g->edge(nc, sub1.second);
    BOOST_REQUIRE(res.second);
    
    LocalEdge et4 = res.first;
    it = g->getGlobalEdges(et4);
    BOOST_CHECK(*it.first == e34);
    BOOST_CHECK(++it.first == it.second);

    //nc to sub3 should have one global edge
    res = g->edge(nc, sub3.second);
    BOOST_REQUIRE(res.second);
    
    LocalEdge et5 = res.first;
    it = g->getGlobalEdges(et5);
    BOOST_CHECK(*it.first == e45);
    BOOST_CHECK(++it.first == it.second);

    //sub1 to sub3 should have one global edge
    res = g->edge(sub1.second, sub3.second);
    BOOST_REQUIRE(res.second);
    
    LocalEdge et6 = res.first;
    it = g->getGlobalEdges(et6);
    BOOST_CHECK(*it.first == e35);
    BOOST_CHECK(++it.first == it.second);

    //v2 to sub2 should have one global edge
    res = g->edge(fusion::at_c<0>(v2), sub2.second);
    BOOST_REQUIRE(res.second);
    
    LocalEdge et7 = res.first;
    it = g->getGlobalEdges(et7);
    BOOST_CHECK(*it.first == e24);
    BOOST_CHECK(++it.first == it.second);
}

BOOST_AUTO_TEST_CASE(filter_graph) {
    
    std::shared_ptr<Graph> g1 = std::shared_ptr<Graph>(new Graph);
 
    fusion::vector<LocalVertex, GlobalVertex> v1c = g1->addVertex();
    fusion::vector<LocalVertex, GlobalVertex> v2c = g1->addVertex();
    fusion::vector<LocalVertex, GlobalVertex> v3c = g1->addVertex();
    fusion::vector<LocalVertex, GlobalVertex> v4c = g1->addVertex();

    fusion::vector<LocalEdge, GlobalEdge, bool> e1c = g1->addEdge(fusion::at_c<0>(v1c),fusion::at_c<0>(v2c));
    fusion::vector<LocalEdge, GlobalEdge, bool> e2c = g1->addEdge(fusion::at_c<0>(v2c),fusion::at_c<0>(v3c));
    fusion::vector<LocalEdge, GlobalEdge, bool> e3c = g1->addEdge(fusion::at_c<0>(v3c),fusion::at_c<0>(v4c));
    
    g1->initIndexMaps();

    dcm::graph::FilterGraph<Graph> filter(g1, 0);
    //at the begining the default edge and vertex group should be 0 and hence the graphs should be identical
    auto eit = g1->edges();
    auto geit = filter.edges();
    for(;eit.first != eit.second; ++eit.first) {
        BOOST_REQUIRE(geit.first != geit.second);
        BOOST_CHECK(*eit.first == *geit.first);
        ++geit.first;
    }
      
    auto vit = g1->vertices();
    auto gvit = filter.vertices();
    for(;vit.first != vit.second; ++vit.first) {
        BOOST_REQUIRE(gvit.first != gvit.second);
        BOOST_CHECK(g1->getProperty<Index>(*vit.first) == filter.getProperty<Index>(*gvit.first));
        ++gvit.first;
    }
   
    //now change two vertices in a new group
    g1->setProperty<Group>(fusion::at_c<0>(v3c), 1);
    g1->setProperty<Group>(fusion::at_c<0>(v4c), 1);
      
    dcm::graph::FilterGraph<Graph> filter2(g1, 1);
 
    auto g2eit = filter2.edges();
    geit = filter.edges();

    //we added two vertices to the new group, so both groups should have two vertices. From the 3 edges we did not
    //shift any. Edges get filtered for their group and the group of the source and target vertices, hence
    //group 1 should ony have one edge remaining (the one between v1 and v2) while group two has non.
    int c = 0;
    for(;geit.first != geit.second; ++geit.first)
        c++;
    BOOST_CHECK(c==1);
    BOOST_CHECK(g2eit.first == g2eit.second);
    
    gvit = filter.vertices();
    auto g2vit = filter2.vertices();
    for(c=0;gvit.first != gvit.second; ++gvit.first)
        c++;
    BOOST_CHECK(c==2);    
    for(c=0;g2vit.first != g2vit.second; ++g2vit.first)
        c++;
    BOOST_CHECK(c==2); 
    
    //now change an edge to group 1
    filter2.setProperty<Group>(fusion::at_c<0>(e3c), 1);
      
    g2eit = filter2.edges();
    geit = filter.edges();

    //we added one edge to the second group. as both vertices of this edge are also in this group
    //it will show up there
    for(c = 0;geit.first != geit.second; ++geit.first)
        c++;
    BOOST_CHECK(c==1);
    for(c=0;g2eit.first != g2eit.second; ++g2eit.first)
        c++;
    BOOST_CHECK(c==1);
    
    gvit = filter.vertices();
    g2vit = filter2.vertices();
    for(c=0;gvit.first != gvit.second; ++gvit.first)
        c++;
    BOOST_CHECK(c==2);    
    for(c=0;g2vit.first != g2vit.second; ++g2vit.first)
        c++;
    BOOST_CHECK(c==2); 
    
}

BOOST_AUTO_TEST_SUITE_END();
