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

#include <opendcm/core/analyse.hpp>

struct value {
    typedef int type;
};

struct symbol {
    typedef int type;
};

typedef dcm::ClusterGraph<mpl::vector1<value>, mpl::vector1<symbol>,
        mpl::vector0<>, mpl::vector0<> > Graph;

struct S1_Node : public dcm::details::TerminalNode<Graph> {
  
    S1_Node(boost::shared_ptr< Graph > g) : dcm::details::TerminalNode<Graph>(g) {};
  
    virtual bool swallow(dcm::GlobalEdge e) {
        if(m_graph->getProperty<symbol>(e) == 1)
            return dcm::details::TerminalNode<Graph>::swallow(e);

        return false;
    };
};

struct S2_Node : public dcm::details::TerminalNode<Graph> {

    S2_Node(boost::shared_ptr< Graph > g) : dcm::details::TerminalNode<Graph>(g) {};
    
    virtual bool swallow(dcm::GlobalEdge e) {
        if(m_graph->getProperty<symbol>(e) == 2)
            return dcm::details::TerminalNode<Graph>::swallow(e);

        return false;
    };
};

BOOST_AUTO_TEST_SUITE(analyse);

BOOST_AUTO_TEST_CASE(analyse_basic) {

    //build up the graph
    boost::shared_ptr< Graph > g = boost::shared_ptr< Graph >(new Graph);
    
    dcm::GlobalVertex v1 = fusion::at_c<1>(g->addVertex());
    dcm::GlobalVertex v2 = fusion::at_c<1>(g->addVertex());
    
    dcm::GlobalEdge e1 = fusion::at_c<1>(g->addEdge(v1,v2));
    dcm::GlobalEdge e2 = fusion::at_c<1>(g->addEdge(v1,v2));
    dcm::GlobalEdge e3 = fusion::at_c<1>(g->addEdge(v1,v2));
    dcm::GlobalEdge e4 = fusion::at_c<1>(g->addEdge(v1,v2));
    
    g->setProperty<symbol>(e1, 1);
    g->setProperty<symbol>(e2, 2);
    g->setProperty<symbol>(e3, 1);
    g->setProperty<symbol>(e4, 2);
    
    //build up the analyser tree
    
    

};

BOOST_AUTO_TEST_SUITE_END();
