/*
    openDCM, dimensional constraint manager
    Copyright (C) 2014  Stefan Troeger <stefantroeger@gmx.net>

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

#define DCM_USE_MULTITHREADING

#include "opendcm/core/scheduler.hpp"
#include <boost/test/unit_test.hpp>
#include <boost/atomic.hpp>
#include <boost/chrono.hpp>

#include <tbb/parallel_for.h>
#include "tbb/blocked_range.h"
#include <tbb/flow_graph.h>

#define NUM_SQRT 1000
#define DCM_DEBUG_THROW_AT_ASSERT

using namespace dcm::details::shedule;

BOOST_AUTO_TEST_SUITE(Scheduler_test_suit);

volatile boost::atomic<int> count;
volatile boost::atomic<int> dummy;
tbb::atomic<int> tbbcount;

struct Test {
       
    void execute() const {  
        volatile int c;
        for ( long i = 0; i < NUM_SQRT; ++i )
           c = std::sqrt(125.34L);// burn some time
            
        count.fetch_add(1, boost::memory_order_relaxed);
    };
};

struct TBBTest {

    void operator() ( int& i ) const {
      
        volatile int c;
        for ( long i = 0; i < NUM_SQRT; ++i )
           c = std::sqrt(125.34L);// burn some time
            
        tbbcount+= 1;
    };
    
    void operator() ( tbb::blocked_range<int>& r ) const {
      
        for ( int i = r.begin(); i != r.end(); ++i )
            operator()(i);
    };
    
    void operator()(const tbb::flow::continue_msg & msg) {
        
        volatile int c;
        for ( long i = 0; i < NUM_SQRT; ++i )
           c = std::sqrt(125.34L);// burn some time
            
        tbbcount+= 1;
    };
};

typedef const tbb::flow::continue_msg & msg_t;
typedef tbb::flow::continue_node< tbb::flow::continue_msg > node_t;

struct TBBGraphBody {

    int nested;
    
    TBBGraphBody(int i) : nested(i) { }
    
    void operator()(const tbb::flow::continue_msg & msg) {
        
        tbb::flow::graph g;
        tbb::flow::broadcast_node<tbb::flow::continue_msg> start(g);
        
        node_t A(g, TBBTest() );
        node_t B(g, TBBTest() );
        node_t C(g, TBBTest() );
        node_t D(g, TBBTest() );
        node_t E(g, TBBTest() );
        node_t F(g, TBBTest() );
        node_t G(g, TBBTest() );
        node_t H(g, TBBTest() );
        node_t I(g, TBBTest() );
        node_t J(g, TBBTest() );
        
        tbb::flow::make_edge(start, A);
        tbb::flow::make_edge(start, B);
        tbb::flow::make_edge(A, C);
        tbb::flow::make_edge(A, D);
        tbb::flow::make_edge(C, E);
        tbb::flow::make_edge(C, F);
        tbb::flow::make_edge(E, G);
        tbb::flow::make_edge(E, H);
        tbb::flow::make_edge(C, I);
        tbb::flow::make_edge(C, J);
        
        if(nested>0) {
            node_t G1(g, TBBGraphBody(nested-1));
            node_t G2(g, TBBGraphBody(nested-1));  
            
            tbb::flow::make_edge(C, G1);
            tbb::flow::make_edge(C, G2);
            
            start.try_put(tbb::flow::continue_msg());
            g.wait_for_all();
        }
        else {
            start.try_put(tbb::flow::continue_msg());
            g.wait_for_all();
        }

        
    };
};
  
Test t;

struct TestGroup : public Group {

    int  nested;
    
    TestGroup(int i) : Group(), nested(i) {
        
        //create some nodes
        std::vector<Node*> nodes(15);
        
        nodes[0] = createNode(boost::bind(&Test::execute, &t));
        nodes[1] = createNode(boost::bind(&Test::execute, &t));
        nodes[2] = createDependendNode(nodes[0], boost::bind(&Test::execute, &t));
        nodes[3] = createDependendNode(nodes[0], boost::bind(&Test::execute, &t));
        nodes[4] = createDependendNode(nodes[2], boost::bind(&Test::execute, &t));
        nodes[5] = createDependendNode(nodes[2], boost::bind(&Test::execute, &t));
        nodes[6] = createDependendNode(nodes[4], boost::bind(&Test::execute, &t));
        nodes[7] = createDependendNode(nodes[4], boost::bind(&Test::execute, &t));
        nodes[8] = createDependendNode(nodes[2], boost::bind(&Test::execute, &t));
        nodes[9] = createDependendNode(nodes[2], boost::bind(&Test::execute, &t));
        
        if(nested>0) {
            TestGroup* group1 = new TestGroup(nested-1);
            TestGroup* group2 = new TestGroup(nested-1);
            setupDependendGroupNode(group1, nodes[2]);
            setupDependendGroupNode(group2, nodes[2]);
        };
    };
};


BOOST_AUTO_TEST_CASE(creation) {
    
    
};

BOOST_AUTO_TEST_CASE(sheduler) {

    typedef boost::chrono::microseconds us;
    typedef boost::chrono::milliseconds ms;
    typedef boost::chrono::seconds s;
    
    try {
        int loop;
        count = 0;
        
        TestGroup* group = new TestGroup(10);
        TBBGraphBody body(10);

        Scheduler sh(2);
        
        boost::chrono::system_clock::time_point start = boost::chrono::system_clock::now();
        sh.execute(group);
        sh.join();
        boost::chrono::duration<double> sec = boost::chrono::system_clock::now() - start;   
           
        std::cout<<count.load(boost::memory_order_relaxed)<<" counts in " << sec.count() << "s"
        << " with 2 threads" << std::endl;
        loop = count;
         
        Scheduler sh2(3);
        count = 0;
        
        start = boost::chrono::system_clock::now();
        sh2.execute(group);
        sh2.join();
        sec = boost::chrono::system_clock::now() - start;   
           
        std::cout<<count.load(boost::memory_order_relaxed)<<" counts in " << sec.count() << "s"
        << " with 3 threads" << std::endl;
        
        Scheduler sh3(4);
        count = 0;
        
        start = boost::chrono::system_clock::now();
        sh3.execute(group);
        sh3.join();
        sec = boost::chrono::system_clock::now() - start;   
           
        std::cout<<count.load(boost::memory_order_relaxed)<<" counts in " << sec.count() << "s"
        << " with 4 threads" << std::endl;
        
        count = 0;
        Test* der = new Test;
        start = boost::chrono::system_clock::now();
        for(int i=0; i<loop; i++)
            der->execute();
        
        sec = boost::chrono::system_clock::now() - start;   
           
        std::cout<<count.load(boost::memory_order_relaxed)<<" counts in " << sec.count() << "s"
        << " by derived" << std::endl;
        
        int c = count.load(boost::memory_order_relaxed);
        tbbcount = 0;
        start = boost::chrono::system_clock::now();
                
        tbb::parallel_for(0,c,1, TBBTest());
        
        sec = boost::chrono::system_clock::now() - start;              
        std::cout<<tbbcount<<" counts in " << sec.count() << "s"
        << " by derived with tbb" << std::endl;
        
        tbbcount = 0;
        start = boost::chrono::system_clock::now();
                
        tbb::parallel_for(tbb::blocked_range<int>(0,c), TBBTest());
        
        sec = boost::chrono::system_clock::now() - start;              
        std::cout<<tbbcount<<" counts in " << sec.count() << "s"
        << " by derived with tbb blocked range" << std::endl;
        
        tbbcount = 0;
        start = boost::chrono::system_clock::now();
                
        body(tbb::flow::continue_msg());
        
        sec = boost::chrono::system_clock::now() - start;              
        std::cout<<tbbcount<<" counts in " << sec.count() << "s"
        << " by tbb dependency" << std::endl;
             
        BOOST_CHECK(count.load(boost::memory_order_relaxed) == 11);

    }
    catch(...) {
        BOOST_FAIL("Exception not expected");
    };
}


BOOST_AUTO_TEST_SUITE_END();
