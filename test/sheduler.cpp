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

#include "opendcm/core/sheduler.hpp"
#include <boost/test/unit_test.hpp>
#include <boost/atomic.hpp>
#include <boost/chrono.hpp>

#define NUM_SQRT 1000
#define DCM_DEBUG_THROW_AT_ASSERT

using namespace dcm::details::shedule;

BOOST_AUTO_TEST_SUITE(Sheduler_test_suit);

volatile boost::atomic<int> count;
volatile boost::atomic<int> dummy;

struct TestNode : public Node {
    bool executed;
    int num;
    
    TestNode(int i) : Node(),executed(false), num(i) {};
    
    virtual void execute() {  
        volatile int c;
        for ( long i = 0; i < NUM_SQRT; ++i )
           c = std::sqrt(125.34L);// burn some time
            
        count.fetch_add(1, boost::memory_order_relaxed);
    };
};

struct TestGroup : public Group {
    bool executed;
    int  nested;
    
    TestGroup(int i) : Group(), executed(false), nested(i) {
        //create some nodes
        std::vector<TestNode*> nodes;
        for(int i=0; i<10; i++)
            nodes.push_back(new TestNode(i));
        
        setupGroupNode(nodes[0]);
        setupGroupNode(nodes[1]);
        setupGroupNode(nodes[2], nodes[0]);
        setupGroupNode(nodes[3], nodes[0]);
        setupGroupNode(nodes[4], nodes[2]);
        setupGroupNode(nodes[5], nodes[2]);
        setupGroupNode(nodes[6], nodes[4]);
        setupGroupNode(nodes[7], nodes[4]);
        setupGroupNode(nodes[8], nodes[2]);
        setupGroupNode(nodes[9], nodes[2]);
        
        if(nested>0) {
            TestGroup* group1 = new TestGroup(nested-1);
            TestGroup* group2 = new TestGroup(nested-1);
            setupGroupNode(group1, nodes[2]);
            setupGroupNode(group2, nodes[2]);
        };
    };
    
    virtual void execute() {
        volatile int c;
        for ( long i = 0; i < NUM_SQRT; ++i )
           c = std::sqrt(125.34L);// burn some time
           
        count.fetch_add(1, boost::memory_order_relaxed);
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
        
        TestGroup* group = new TestGroup(6);

        Sheduler sh(1);
        
        boost::chrono::system_clock::time_point start = boost::chrono::system_clock::now();
        sh.execute(group);
        sh.join();
        boost::chrono::duration<double> sec = boost::chrono::system_clock::now() - start;   
           
        std::cout<<count.load(boost::memory_order_relaxed)<<" counts in " << sec.count() << "s"
        << " with 1 thread" << std::endl;
        loop = count;
         
        Sheduler sh2(2);
        count = 0;
        
        start = boost::chrono::system_clock::now();
        sh2.execute(group);
        sh2.join();
        sec = boost::chrono::system_clock::now() - start;   
           
        std::cout<<count.load(boost::memory_order_relaxed)<<" counts in " << sec.count() << "s"
        << " with 2 threads" << std::endl;
        
        Sheduler sh3(3);
        count = 0;
        
        start = boost::chrono::system_clock::now();
        sh3.execute(group);
        sh3.join();
        sec = boost::chrono::system_clock::now() - start;   
           
        std::cout<<count.load(boost::memory_order_relaxed)<<" counts in " << sec.count() << "s"
        << " with 3 threads" << std::endl;
        
        count = 0;
        Node* der = new TestNode(1);
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
