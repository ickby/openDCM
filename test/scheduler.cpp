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

//#define DCM_USE_MULTITHREADING

#include "opendcm/core/scheduler.hpp"
#include <boost/test/unit_test.hpp>
#include <boost/atomic.hpp>
#include <boost/chrono.hpp>

#define NUM_SQRT 2000
#define DCM_DEBUG_THROW_AT_ASSERT

using namespace dcm::details::shedule;

BOOST_AUTO_TEST_SUITE(Scheduler_test_suit);

volatile boost::atomic<int> count;
volatile boost::atomic<int> dummy;

struct Test {
       
    void execute() {  
        volatile int c;
        for ( long i = 0; i < NUM_SQRT; ++i )
           c = std::sqrt(125.34L);// burn some time
            
        count.fetch_add(1, boost::memory_order_relaxed);
    };
};

Test t;

struct TestGroup : public Group {

    int  nested;
    
    TestGroup(int i) : Group(), nested(i) {
        
        //create some nodes
        std::vector<Node*> nodes(10);
        
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
        
        TestGroup* group = new TestGroup(8);

        Scheduler sh(1);
        
        boost::chrono::system_clock::time_point start = boost::chrono::system_clock::now();
        sh.execute(group);
        sh.join();
        boost::chrono::duration<double> sec = boost::chrono::system_clock::now() - start;   
           
        std::cout<<count.load(boost::memory_order_relaxed)<<" counts in " << sec.count() << "s"
        << " with 1 thread" << std::endl;
        loop = count;
         
        Scheduler sh2(2);
        count = 0;
        
        start = boost::chrono::system_clock::now();
        sh2.execute(group);
        sh2.join();
        sec = boost::chrono::system_clock::now() - start;   
           
        std::cout<<count.load(boost::memory_order_relaxed)<<" counts in " << sec.count() << "s"
        << " with 2 threads" << std::endl;
        
        Scheduler sh3(3);
        count = 0;
        
        start = boost::chrono::system_clock::now();
        sh3.execute(group);
        sh3.join();
        sec = boost::chrono::system_clock::now() - start;   
           
        std::cout<<count.load(boost::memory_order_relaxed)<<" counts in " << sec.count() << "s"
        << " with 3 threads" << std::endl;
        
        count = 0;
        Test* der = new Test;
        start = boost::chrono::system_clock::now();
        for(int i=0; i<loop; i++)
            der->execute();
        
        sec = boost::chrono::system_clock::now() - start;   
           
        std::cout<<count.load(boost::memory_order_relaxed)<<" counts in " << sec.count() << "s"
        << " by derived" << std::endl;
             
        BOOST_CHECK(count.load(boost::memory_order_relaxed) == 11);

    }
    catch(...) {
        BOOST_FAIL("Exception not expected");
    };
}


BOOST_AUTO_TEST_SUITE_END();
