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
/*
#include "opendcm/core/object.hpp"

#include <boost/function.hpp>
#include <boost/bind.hpp>

#include <boost/mpl/map.hpp>

#ifdef DCM_EXTERNAL_CORE
#include "opendcm/core/imp/system_imp.hpp"
#include "opendcm/core/imp/clustergraph_imp.hpp"
#include "opendcm/core/imp/object_imp.hpp"
#endif
*/
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(system_and_object);
/*
struct test_edge_property1 {
    typedef dcm::edge_property kind;
    typedef int type;
};
struct test_vertex_property1 {
    typedef dcm::vertex_property kind;
    typedef int type;
};
struct test_edge_property2 {
    typedef dcm::edge_property kind;
    typedef int type;
};
struct test_vertex_property2 {
    typedef dcm::vertex_property kind;
    typedef int type;
};
struct test_signal1 {};
struct test_signal2 {};
struct test_signal3 {};
*/
struct TestModule1 {

    typedef boost::mpl::int_<1> ID;

    template<typename Final, typename Stacked>
    struct type : public Stacked {

	int module_function1() {return 1;};
      
        struct TestType1 {
            int value1;
            int function1() {
                return 1;
            };
        };
    };
};
struct TestModule2 {

    typedef boost::mpl::int_<2> ID;

    template<typename Final, typename Stacked>
    struct type : public Stacked {

	int module_function2() {return 2;};
      
        struct TestType1 : public Stacked::TestType1 {
            int value2;
            int function2() {
                return 2;
            };
        };
    };
};

typedef dcm::System<TestModule2, TestModule1> System;



BOOST_AUTO_TEST_CASE(inherit_functions_types) {

    System sys;
    
    BOOST_CHECK(sys.module_function1() == 1);
    BOOST_CHECK(sys.module_function2() == 2);

};
/*
BOOST_AUTO_TEST_CASE(graph_properties) {

    System sys;

    dcm::GlobalVertex v = fusion::at_c<1>(sys.m_cluster->addVertex());
    sys.m_cluster->setProperty<test_vertex_property1>(v, 1);
    sys.m_cluster->setProperty<test_vertex_property2>(v, 2);
    BOOST_CHECK(sys.m_cluster->getProperty<test_vertex_property1>(v) == 1);
    BOOST_CHECK(sys.m_cluster->getProperty<test_vertex_property2>(v) == 2);

    dcm::GlobalVertex v2 = fusion::at_c<1>(sys.m_cluster->addVertex());
    dcm::GlobalEdge e = fusion::at_c<1>(sys.m_cluster->addEdge(v, v2));
    sys.m_cluster->setProperty<test_edge_property1>(e, 1);
    sys.m_cluster->setProperty<test_edge_property2>(e, 2);
    BOOST_CHECK(sys.m_cluster->getProperty<test_edge_property1>(e) == 1);
    BOOST_CHECK(sys.m_cluster->getProperty<test_edge_property2>(e) == 2);

};

BOOST_AUTO_TEST_CASE(object_properties) {

    System sys;
    typedef TestModule1::type<System> Module1;
    typedef TestModule2::type<System> Module2;
    typedef Module1::test_object1 to1;
    typedef Module2::test_object2 to2;

    to1 o1(sys);
    to2 o2(sys);

    //check defualt value
    BOOST_CHECK(o1.getProperty<Module1::test_object1_prop>() == 3);

    o1.setProperty<Module1::test_object1_prop>(5);
    BOOST_CHECK(o1.getProperty<Module1::test_object1_prop>() == 5);

    o1.setProperty<Module2::test_object1_external_prop>(7);
    BOOST_CHECK(o1.getProperty<Module1::test_object1_prop>() == 5);
    BOOST_CHECK(o1.getProperty<Module2::test_object1_external_prop>() == 7);

    o2.getProperty<Module2::test_object2_prop>() = o1.getProperty<Module2::test_object1_external_prop>();
    BOOST_CHECK(o2.getProperty<Module2::test_object2_prop>() == 7);

};

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

