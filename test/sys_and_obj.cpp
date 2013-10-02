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
#include "opendcm/core/object.hpp"

#include <boost/function.hpp>
#include <boost/bind.hpp>

#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(system_and_object);

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

struct TestModule1 {

    template<typename Sys>
    struct type {
        typedef mpl::map< mpl::pair<test_signal1, boost::function<void ()> >,
                mpl::pair<test_signal2, boost::function<void (double, double)> > > signal_map;

        typedef dcm::Unspecified_Identifier Identifier;

        template<typename Syst, typename Derived>
        struct test_object1_base : public dcm::Object<Syst, Derived, signal_map > {
            test_object1_base(Syst& system) : dcm::Object<Syst, Derived, signal_map >(system) { };
            int value;

            void emit_test_void() {
                dcm::Object<Syst, test_object1, signal_map >::template emitSignal<test_signal1>();
            };
            void emit_test_double(const double& d1, const double& d2) {
                dcm::Object<Syst, test_object1, signal_map >::template emitSignal<test_signal2>(d1, d2);
            };
        };

        struct test_object1 : public test_object1_base<Sys, test_object1> {
            test_object1(Sys& system) : test_object1_base<Sys, test_object1>(system) { };
        };

        struct inheriter {
            int test_inherit1() {
                return 1;
            };
            int test_inherit2(int f, int s) {
                return f+s;
            };
        };

        struct test_object1_prop {
            typedef int type;
            typedef test_object1 kind;
            struct default_value {
                int operator()() {
                    return 3;
                };
            };
        };
	
	struct setting1_prop {
            typedef bool type;
            typedef dcm::setting_property kind;
            struct default_value {
                int operator()() {
                    return true;
                };
            };
        };

        typedef mpl::vector0<> geometries;
        typedef mpl::vector1<test_object1> objects;
        typedef mpl::vector4<test_edge_property1, test_vertex_property1, test_object1_prop, setting1_prop>   properties;

        template<typename System>
        static void system_init(System& sys) {};
    };
};

struct TestModule2 {

    template<typename Sys>
    struct type {
        typedef mpl::map< mpl::pair<test_signal1, boost::function<void ()> >,
                mpl::pair<test_signal2, boost::function<void (double, double)> > > signal_map;

        struct test_object2 : public dcm::Object<Sys, test_object2, signal_map > {
            test_object2(Sys& system) : dcm::Object<Sys, test_object2, signal_map >(system) { };
            int value;
        };

        struct inheriter {
            template<typename T1>
            int test_inherit3(T1 f, T1 s) {
                return ((Sys*)this)->test_inherit2(f,s);
            };
        };

        struct test_object1_external_prop {
            typedef int type;
            typedef typename TestModule1::type<Sys>::test_object1 kind;
        };
        struct test_object2_prop {
            typedef int type;
            typedef test_object2 kind;
        };

        typedef mpl::vector1<test_object2>  	objects;
        typedef mpl::vector0<> geometries;
        typedef mpl::vector4<test_edge_property2, test_vertex_property2,
                test_object2_prop, test_object1_external_prop> properties;
        typedef dcm::Unspecified_Identifier Identifier;

        template<typename System>
        static void system_init(System& sys) {};
    };
};

typedef dcm::System<dcm::Kernel<double>, TestModule1, TestModule2> System;

BOOST_AUTO_TEST_CASE(inherit_functions) {

    System sys;
    BOOST_CHECK(sys.test_inherit1() == 1);
    BOOST_CHECK(sys.test_inherit2(2,3) == 5);
    BOOST_CHECK(sys.test_inherit3(2,3) == 5);
};

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
    BOOST_CHECK(sys.getSetting<Module1::setting1_prop>());
    sys.setSetting<Module1::setting1_prop>(false);
    BOOST_CHECK(!sys.getSetting<Module1::setting1_prop>());
    sys.getSetting<Module1::setting1_prop>() = true;
    BOOST_CHECK(sys.getSetting<Module1::setting1_prop>());
    
    //test kernel settings
    BOOST_CHECK(sys.setting<dcm::precision>() == 1e-6);
    sys.setting<dcm::precision>() = 4e-7;
    BOOST_CHECK(sys.setting<dcm::precision>() == 4e-7);
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

BOOST_AUTO_TEST_CASE(object_signals) {

    System sys;

    typedef TestModule1::type<System> Module1;
    typedef Module1::test_object1 to1;

    to1 o1(sys);

    test_functor_void s, s2;

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
};

BOOST_AUTO_TEST_SUITE_END();

