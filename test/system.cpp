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

#include "opendcm/core.hpp"


#include <boost/function.hpp>
#include <boost/bind.hpp>

#include <boost/mpl/map.hpp>
#include <boost/mpl/replace_if.hpp>

#include <boost/test/unit_test.hpp>
#include <boost/exception/get_error_info.hpp>

BOOST_AUTO_TEST_SUITE(system_and_object);

struct TestProperty1 {
    typedef int type;
    struct default_value {
        int operator()() {
            return 3;
        }
    };
};
struct TestProperty2 {
    typedef int type;
    struct change_tracking {};
};
struct TestProperty3 {
    typedef int type;
};
struct TestProperty4 {
    typedef bool type;
};
struct TestVertexProperty1 {
    typedef int type;
};
struct TestVertexProperty2 {
    typedef int type;
};
struct TestEdgeProperty1 {
    typedef int type;
};
struct TestEdgeProperty2 {
    typedef int type;
};
struct TestEdgeProperty3 {
    typedef int type;
};


struct test_signal1 {};
struct test_signal2 {};
struct test_signal3 {};

struct TestModule1 {

    typedef boost::mpl::int_<1> ID;

    template<typename Final, typename Stacked>
    struct type : public Stacked {

        int module_function1() {
            return 1;
        };

        struct TestType1 : public Stacked::Object {

            TestType1() : Stacked::Object( Final::template objectTypeID<typename Final::TestType1>::ID::value ) {};
            int function1() {
                return 1;
            };
            
            typedef typename Stacked::Object Base;
            DCM_OBJECT_ADD_PROPERTIES( Base, (TestProperty1)(TestProperty2) )
        };
        
        struct TestType2 : public Stacked::Object {
                    
            TestType2() : Stacked::Object( Final::template objectTypeID<typename Final::TestType2>::ID::value ) {};
                        
            typedef typename Stacked::Object Base;
            DCM_OBJECT_ADD_PROPERTIES( Base, (TestProperty3) )
        };

        DCM_MODULE_ADD_OBJECTS(Stacked, (TestType1)(TestType2))
        DCM_MODULE_ADD_VERTEX_PROPERTIES(Stacked, (TestVertexProperty1))
        DCM_MODULE_ADD_GLOBAL_EDGE_PROPERTIES(Stacked, (TestEdgeProperty3))
    };
};

struct TestModule2 {

    typedef boost::mpl::int_<2> ID;

    template<typename Final, typename Stacked>
    struct type : public Stacked {

        int module_function2() {
            return 2;
        };

        struct TestType1 : public Stacked::TestType1 {
            int function2() {
                return 2;
            };
            
            typedef typename Stacked::TestType1 Base;
            DCM_OBJECT_ADD_PROPERTIES(Base, (TestProperty4))
        };
        
        DCM_MODULE_ADD_OBJECTS(Stacked, (TestType1))
        DCM_MODULE_ADD_VERTEX_PROPERTIES(Stacked, (TestVertexProperty2))
        DCM_MODULE_ADD_LOCAL_EDGE_PROPERTIES(Stacked, (TestEdgeProperty1)(TestEdgeProperty2))
    };
};

typedef dcm::System<TestModule2, TestModule1> System;
typedef dcm::System<dcm::Eigen3Kernel<double>, TestModule2, TestModule1> System2;



BOOST_AUTO_TEST_CASE(module_functions) {

    System2 sys;

    BOOST_CHECK(sys.module_function1() == 1);
    BOOST_CHECK(sys.module_function2() == 2);

};

BOOST_AUTO_TEST_CASE(object_handling) {

    try {
        
    boost::shared_ptr<System::TestType1>  t(new System::TestType1);
    boost::shared_ptr<System::Object> b(t);
    
    //test the type ID stuff
    BOOST_CHECK(t->isType<System::TestType1>());
    BOOST_CHECK(!t->isType<System::TestType2>());
    
    //test the property stuff
    BOOST_CHECK(b->isType<System::TestType1>());
    BOOST_CHECK(!b->isType<System::TestType2>());

    BOOST_CHECK(t->getProperty<TestProperty1>() == 3);
    t->setProperty<TestProperty1>(2);
    BOOST_CHECK(t->getProperty<TestProperty1>() == 2);
    t->setProperty<TestProperty4>(true);
    BOOST_CHECK(t->getProperty<TestProperty4>());

    /*
    b->setProperty<TestProperty1>(5);
    BOOST_CHECK(b->getProperty<TestProperty1>() == 5);
    b->setProperty<TestProperty4>(false);
    BOOST_CHECK(!b->getProperty<TestProperty4>());*/
    
    t->setTestProperty1(2);
    BOOST_CHECK(t->getTestProperty1());
    t->setTestProperty4(true);
    BOOST_CHECK(t->getTestProperty4());
   
    BOOST_CHECK(t->function1() == 1);
    BOOST_CHECK(t->function2() == 2);
    
    }
    catch(boost::exception& x) {
        BOOST_FAIL(*boost::get_error_info<dcm::error_message>(x));
    }
    catch(std::exception& x) {
        BOOST_FAIL("Unknown exception");
    }
};

BOOST_AUTO_TEST_CASE(graph_properties) {

    System sys;
    auto gr = std::static_pointer_cast<typename System::Graph>(sys.getGraph());
    
    dcm::graph::GlobalVertex v = fusion::at_c<1>(gr->addVertex());
    gr->setProperty<TestVertexProperty1>(v, 1);
    gr->setProperty<TestVertexProperty2>(v, 2);
    BOOST_CHECK(gr->getProperty<TestVertexProperty1>(v) == 1);
    BOOST_CHECK(gr->getProperty<TestVertexProperty2>(v) == 2);

    dcm::graph::GlobalVertex v2 = fusion::at_c<1>(gr->addVertex());
    auto ret = gr->addEdge(v, v2);
    dcm::graph::GlobalEdge e = fusion::at_c<1>(ret);
    dcm::graph::LocalEdge le = fusion::at_c<0>(ret);
    gr->setProperty<TestEdgeProperty3>(e, 1);
    gr->setProperty<TestEdgeProperty1>(le, 2);
    gr->setProperty<TestEdgeProperty2>(le, 3);
    BOOST_CHECK(gr->getProperty<TestEdgeProperty3>(e) == 1);
    BOOST_CHECK(gr->getProperty<TestEdgeProperty1>(le) == 2);
    BOOST_CHECK(gr->getProperty<TestEdgeProperty2>(le) == 3);

};

/*
BOOST_AUTO_TEST_CASE(settings_properties) {

    System sys;
    typedef TestModule1::type<System> Module1;

    //check general settings
    BOOST_CHECK(sys.getOption<Module1::setting1_prop>());
    sys.setOption<Module1::setting1_prop>(false);
    BOOST_CHECK(!sys.getOption<Module1::setting1_prop>());
    sys.getOption<Module1::setting1_prop>() = true;
    BOOST_CHECK(sys.getOption<Module1::setting1_prop>());

    //test kernel settings
    BOOST_CHECK(sys.option<dcm::precision>() == 1e-6);
    sys.option<dcm::precision>() = 4e-7;
    BOOST_CHECK(sys.option<dcm::precision>() == 4e-7);
};

struct test_functor_void {
    test_functor_void() : counter(0) {};
    void count() {
        counter++;
    };
    int counter;
};

struct test_functor_double {
    test_functor_double() : counter(0) {};
    void count(double d1, double d2) {
        counter += int(d1+d2);
    };
    int counter;
};

BOOST_AUTO_TEST_CASE(signals) {

    System sys;

    typedef TestModule1::type<System> Module1;
    typedef Module1::test_object1 to1;

    to1 o1(sys);

    test_functor_void s, s2, s3;

    dcm::Connection c1 = o1.connectSignal<test_signal1>(boost::bind(&test_functor_void::count, &s));
    o1.connectSignal<test_signal1>(boost::bind(&test_functor_void::count, &s2));

    o1.emit_test_void();
    o1.emit_test_void();

    BOOST_CHECK(s.counter == 2);
    BOOST_CHECK(s2.counter == 2);

    o1.disconnectSignal<test_signal1>(c1);
    o1.emit_test_void();

    BOOST_CHECK(s.counter == 2);
    BOOST_CHECK(s2.counter == 3);

    test_functor_double d;

    o1.connectSignal<test_signal2>(boost::bind(&test_functor_double::count, &d, _1, _2));
    o1.emit_test_double(2,4);

    BOOST_CHECK(d.counter == 6);

    sys.connectSignal<test_signal3>(boost::bind(&test_functor_void::count, &s3));
    sys.emit_signal_3();
    sys.emit_signal_3();

    BOOST_CHECK(s3.counter == 2);
};*/

BOOST_AUTO_TEST_SUITE_END();

