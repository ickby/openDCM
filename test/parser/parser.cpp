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

#include "parser.hpp"

#include "opendcm/moduleState/indent.hpp"
#include <boost/iostreams/filtering_stream.hpp>

//we need isSame implementation
#include "opendcm/core/imp/kernel_imp.hpp"

template<typename T, typename G>
void generate(std::stringstream& s, T& input, G& gen) {

    boost::iostreams::filtering_ostream indent_stream;
    indent_stream.push(indent_filter());
    indent_stream.push(s);

    std::ostream_iterator<char> out(indent_stream);
    karma::generate(out, gen, input);
};
template<typename T, typename G>
void parse(std::istream& stream, T& output, G& par) {

    //disable skipping of whitespace
    stream.unsetf(std::ios::skipws);

    // wrap istream into iterator
    boost::spirit::istream_iterator begin(stream);
    boost::spirit::istream_iterator end;

    // use iterator to parse file data
    qi::phrase_parse(begin, end, par, qi::space, output);
};

BOOST_AUTO_TEST_SUITE(parser_suit);
/*
BOOST_AUTO_TEST_CASE(parser_seperate) {

    //test properties
    std::stringstream s1;
    dcm::details::prop_grammar<dcm::type_prop,
        dcm::parser_generator<dcm::type_prop, System, std::ostream_iterator<char> > > gram;
    int value = 5;
    generate(s1, value, gram);


    int result;
    dcm::details::prop_parser<dcm::type_prop,
        dcm::parser_parser<dcm::type_prop, System, boost::spirit::istream_iterator > > pars;
    parse(s1, result, pars);

    BOOST_CHECK(result == 5);

    //test property vector
    std::stringstream s2;
    dcm::details::cluster_prop_gen<System> gram2;
    System::Cluster::Properties bundle;
    fusion::at_c<0>(bundle) = false;
    fusion::at_c<1>(bundle) = 15;
    generate(s2, bundle, gram2);

    System::Cluster::Properties bundle2;
    dcm::details::cluster_prop_par<System> pars2;
    parse(s2, bundle2, pars2);

    BOOST_CHECK(fusion::at_c<0>(bundle2) == false);
    BOOST_CHECK(fusion::at_c<1>(bundle2) == 15);
    //std::cout<<s2.str()<<std::endl;
    //std::cout<<"result: " << result<<std::endl;
};
*/
BOOST_AUTO_TEST_CASE(parser_graph) {

    System sys;

    fusion::vector<dcm::LocalVertex, dcm::GlobalVertex> res1 = sys.m_cluster->addVertex();
    fusion::vector<dcm::LocalVertex, dcm::GlobalVertex> res2 = sys.m_cluster->addVertex();
    fusion::vector<dcm::LocalVertex, dcm::GlobalVertex> res3 = sys.m_cluster->addVertex();
    dcm::LocalEdge e1 = fusion::at_c<0>(sys.m_cluster->addEdge(fusion::at_c<1>(res1),fusion::at_c<1>(res2)));
    dcm::LocalEdge e2 = fusion::at_c<0>(sys.m_cluster->addEdge(fusion::at_c<1>(res2),fusion::at_c<1>(res3)));
    dcm::LocalEdge e3 = fusion::at_c<0>(sys.m_cluster->addEdge(fusion::at_c<1>(res1),fusion::at_c<1>(res3)));

    boost::shared_ptr<TestModule1::type<System>::test_object1> ptr(new TestModule1::type<System>::test_object1(sys));
    sys.m_cluster->setObject(fusion::at_c<1>(res1), ptr);

    sys.m_cluster->setProperty<TestModule1::type<System>::test_vertex1_prop>(fusion::at_c<0>(res1), 1);
    sys.m_cluster->setProperty<TestModule1::type<System>::test_vertex1_prop>(fusion::at_c<0>(res2), 2);
    sys.m_cluster->setProperty<TestModule1::type<System>::test_vertex1_prop>(fusion::at_c<0>(res3), 3);
    sys.m_cluster->setProperty<TestModule1::type<System>::test_edge1_prop>(e1, 1);
    sys.m_cluster->setProperty<TestModule1::type<System>::test_edge1_prop>(e2, 2);
    sys.m_cluster->setProperty<TestModule1::type<System>::test_edge1_prop>(e3, 3);

    //subcluster
    boost::shared_ptr<System::Cluster> scl1 = sys.m_cluster->createCluster().first;
    boost::shared_ptr<System::Cluster> scl2 = sys.m_cluster->createCluster().first;

    fusion::vector<dcm::LocalVertex, dcm::GlobalVertex> res4 = scl1->addVertex();
    fusion::vector<dcm::LocalVertex, dcm::GlobalVertex> res5 = scl1->addVertex();
    fusion::vector<dcm::LocalVertex, dcm::GlobalVertex> res6 = scl2->addVertex();
    dcm::LocalEdge e4 = fusion::at_c<0>(sys.m_cluster->addEdge(fusion::at_c<1>(res4),fusion::at_c<1>(res5)));
    dcm::LocalEdge e5 = fusion::at_c<0>(sys.m_cluster->addEdge(fusion::at_c<1>(res5),fusion::at_c<1>(res6)));
    dcm::LocalEdge e6 = fusion::at_c<0>(sys.m_cluster->addEdge(fusion::at_c<1>(res1),fusion::at_c<1>(res4)));

    scl1->setProperty<TestModule1::type<System>::test_vertex1_prop>(fusion::at_c<0>(res4), 4);
    scl1->setProperty<TestModule1::type<System>::test_vertex1_prop>(fusion::at_c<0>(res5), 5);
    scl2->setProperty<TestModule1::type<System>::test_vertex1_prop>(fusion::at_c<0>(res6), 6);
    scl1->setProperty<TestModule1::type<System>::test_edge1_prop>(e4, 7);

    sys.m_cluster->setProperty<dcm::type_prop>(5);
    sys.m_cluster->setProperty<dcm::changed_prop>(false);
    scl1->setProperty<dcm::type_prop>(2);
    scl1->setProperty<dcm::changed_prop>(false);
    scl2->setProperty<dcm::type_prop>(7);
    scl2->setProperty<dcm::changed_prop>(true);

    std::stringstream s;
    sys.saveState(s);

    //std::cout<<s.str()<<std::endl;
     
    //change main graphs property as this is not reset
    sys.m_cluster->setProperty<dcm::type_prop>(2);
    sys.m_cluster->setProperty<dcm::changed_prop>(false);

    sys.clear();
    BOOST_CHECK(boost::num_vertices(*sys.m_cluster) == 0);
    BOOST_CHECK(boost::num_edges(*sys.m_cluster) == 0);
    BOOST_REQUIRE(sys.m_cluster->numClusters() == 0);
    
    //change properties so we can see if they get rewritten
    sys.m_cluster->getProperty<dcm::type_prop>() = 35;
    sys.m_cluster->getProperty<dcm::changed_prop>() = true;
    scl1->getProperty<dcm::type_prop>() = 32;
    scl2->getProperty<dcm::type_prop>() = 37;
      
    //load the state
    sys.loadState(s);
 
    //check clusters
    BOOST_CHECK(boost::num_vertices(*sys.m_cluster) == 5);
    BOOST_CHECK(boost::num_edges(*sys.m_cluster) == 5);
    BOOST_REQUIRE(sys.m_cluster->numClusters() == 2);
    
    //access the two subcluster
    typedef System::Cluster::cluster_iterator cit;
    std::pair<cit, cit> it = sys.m_cluster->clusters();
    boost::shared_ptr<System::Cluster> nscl1 = (*(it.first)).second;
    boost::shared_ptr<System::Cluster> nscl2 = (*(++it.first)).second;
    
    BOOST_CHECK(sys.m_cluster->getProperty<dcm::type_prop>() == 5);
    BOOST_CHECK(!sys.m_cluster->getProperty<dcm::changed_prop>());
    BOOST_CHECK(nscl1->getProperty<dcm::type_prop>() == 2);
    BOOST_CHECK(!nscl1->getProperty<dcm::changed_prop>());
    BOOST_CHECK(nscl2->getProperty<dcm::type_prop>() == 7);
    BOOST_CHECK(nscl2->getProperty<dcm::changed_prop>());

    //check vertex and edge properties
    std::pair<dcm::LocalVertex, bool> v1 = sys.m_cluster->getLocalVertex(fusion::at_c<1>(res1));
    std::pair<dcm::LocalVertex, bool> v2 = sys.m_cluster->getLocalVertex(fusion::at_c<1>(res2));
    std::pair<dcm::LocalVertex, bool> v3 = sys.m_cluster->getLocalVertex(fusion::at_c<1>(res3));

    BOOST_REQUIRE(v1.second);
    BOOST_REQUIRE(v2.second);
    BOOST_REQUIRE(v3.second);

    BOOST_CHECK(sys.m_cluster->getProperty<TestModule1::type<System>::test_vertex1_prop>(v1.first) == 1);
    BOOST_CHECK(sys.m_cluster->getProperty<TestModule1::type<System>::test_vertex1_prop>(v2.first) == 2);
    BOOST_CHECK(sys.m_cluster->getProperty<TestModule1::type<System>::test_vertex1_prop>(v3.first) == 3);
    BOOST_CHECK(sys.m_cluster->getProperty<TestModule1::type<System>::test_edge1_prop>(sys.m_cluster->edge(v1.first, v2.first).first) == 1);
    BOOST_CHECK(sys.m_cluster->getProperty<TestModule1::type<System>::test_edge1_prop>(sys.m_cluster->edge(v2.first,v3.first).first) == 2);
    BOOST_CHECK(sys.m_cluster->getProperty<TestModule1::type<System>::test_edge1_prop>(sys.m_cluster->edge(v1.first,v3.first).first) == 3);

    std::pair<dcm::LocalVertex, bool> v4 = nscl1->getLocalVertex(fusion::at_c<1>(res4));
    std::pair<dcm::LocalVertex, bool> v5 = nscl1->getLocalVertex(fusion::at_c<1>(res5));
    std::pair<dcm::LocalVertex, bool> v6 = nscl2->getLocalVertex(fusion::at_c<1>(res6));

    BOOST_REQUIRE(v4.second);
    BOOST_REQUIRE(v5.second);
    BOOST_REQUIRE(v6.second);

    BOOST_CHECK(nscl1->getProperty<TestModule1::type<System>::test_vertex1_prop>(v4.first) == 4);
    BOOST_CHECK(nscl1->getProperty<TestModule1::type<System>::test_vertex1_prop>(v5.first) == 5);
    BOOST_CHECK(nscl2->getProperty<TestModule1::type<System>::test_vertex1_prop>(v6.first) == 6);
    BOOST_CHECK(nscl1->getProperty<TestModule1::type<System>::test_edge1_prop>(nscl1->edge(v4.first, v5.first).first) == 7);

}

BOOST_AUTO_TEST_CASE(parser_module3d) {

  System sys, sys2;
  
  Eigen::Vector3d v1, v2;
  v1 << 1,2,3;
  v2 << 1.23456789, 2.34567891, 3.45678912;
  
  boost::shared_ptr<Geometry3D> g1 = sys.createGeometry3D(v1, 1);
  boost::shared_ptr<Geometry3D> g2 = sys.createGeometry3D(v2, 2);
  sys.createConstraint3D(3, g1, g2, dcm::distance=3.);
  
  sys.setOption<dcm::precision>(2e-9);
  
  std::stringstream s;
  sys.saveState(s);
  
  //some things are note reset, to test we need to override them first!
  sys.setOption<dcm::precision>(5e-2);
  
  std::cout<<s.str()<<std::endl;
  
  sys2.loadState(s);
  
  BOOST_REQUIRE( sys2.hasGeometry3D(1) );
  BOOST_REQUIRE( sys2.hasGeometry3D(2) );
  BOOST_REQUIRE( sys2.hasConstraint3D(3) );
  
  boost::shared_ptr<Geometry3D> ng1 = sys2.getGeometry3D(1);
  boost::shared_ptr<Geometry3D> ng2 = sys2.getGeometry3D(2);
  boost::shared_ptr<Constraint3D> nc1 = sys2.getConstraint3D(3);
   
  BOOST_REQUIRE( ng1 );
  BOOST_REQUIRE( ng2 );
  BOOST_REQUIRE( nc1 );
   
  BOOST_CHECK(nc1->first == ng1);
  BOOST_CHECK(nc1->second== ng2);

  Eigen::Vector3d nv1, nv2;
  nv1 = get<Eigen::Vector3d>(ng1);
  nv2 = get<Eigen::Vector3d>(ng2);
  
  BOOST_CHECK( nv1.isApprox(v1) );
  BOOST_CHECK( nv2.isApprox(v2) );
  
  sys2.solve();
  
  nv1 = get<Eigen::Vector3d>(ng1);
  nv2 = get<Eigen::Vector3d>(ng2);
  
  BOOST_CHECK( Kernel::isSame((nv1-nv2).norm(), 3, 10e-6) );
  BOOST_CHECK( Kernel::isSame(sys.getOption<dcm::precision>(), 2e-9, 1e-12) );
  
  std::cout<<"precision: "<<sys.getOption<dcm::precision>()<<std::endl;
  
}

BOOST_AUTO_TEST_SUITE_END();
