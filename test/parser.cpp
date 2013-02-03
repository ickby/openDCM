/*
    openDCM, dimensional constraint manager
    Copyright (C) 2013  Stefan Troeger <stefantroeger@gmx.net>

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

#include "opendcm/Core"
#include "opendcm/ModuleState"
#include "opendcm/core/property.hpp"

#include <iosfwd>
#include <sstream>
#include <boost/spirit/include/karma_int.hpp>

#include <boost/test/unit_test.hpp>

using boost::spirit::karma::lit;
using boost::spirit::ascii::string;
using boost::spirit::karma::int_;


struct TestModule1 {

    template<typename Sys>
    struct type {
        typedef mpl::map0<> signal_map;

        struct test_object1 : public dcm::Object<Sys, test_object1, signal_map > {
            test_object1(Sys& system) : dcm::Object<Sys, test_object1, signal_map >(system) { };
        };
        struct test_object2 : public dcm::Object<Sys, test_object2, signal_map > {
            test_object2(Sys& system) : dcm::Object<Sys, test_object2, signal_map >(system) { };
        };

        struct inheriter {};

        struct test_object1_prop {
            typedef int type;
            typedef test_object1 kind;
        };
        struct test_object2_prop {
            typedef std::string type;
            typedef test_object1 kind;
        };
        struct test_vertex1_prop {
            typedef std::string type;
            typedef dcm::vertex_property kind;
        };
        struct test_edge1_prop {
            typedef int type;
            typedef dcm::edge_property kind;
        };

        typedef mpl::vector2<test_object1, test_object2> objects;
        typedef mpl::vector4<test_object1_prop, test_object2_prop, test_vertex1_prop, test_edge1_prop>   properties;

        static void system_init(Sys& sys) {};
    };
};

struct test_prop {
    typedef std::string type;
    typedef int kind;
};

namespace dcm {
/*
template<typename System, typename iterator>
struct parser_generator< typename TestModule1::type<System>::test_object1_prop, System, iterator > {

    typedef rule<iterator, int()> generator;
    static void init(generator& r) {
        r = lit("<type>test o1 property</type>\n<value>")<<int_<<"</value>";
    };
};
*/
template<typename System, typename iterator>
struct parser_generator< test_prop, System, iterator > {

    typedef rule<iterator, std::string()> generator;
    static void init(generator& r) {
        r = lit("<type>test property</type>\n<value>")<<string<<"</value>";
    };
};
template<typename System, typename iterator>
struct parser_generator< typename TestModule1::type<System>::test_vertex1_prop, System, iterator > {

    typedef rule<iterator, std::string()> generator;
    static void init(generator& r) {
        r = lit("<type>vertex 1 prop</type>\n<value>")<<string<<"</value>";
    };
};
template<typename System, typename iterator>
struct parser_generator< typename TestModule1::type<System>::test_edge1_prop, System, iterator > {

    typedef rule<iterator, int()> generator;
    static void init(generator& r) {
        r = lit("<type>vedge1 prop</type>\n<value>")<<int_<<"</value>";
    };
};

template<typename System>
struct parser_generate< typename TestModule1::type<System>::test_object1, System>
  : public mpl::true_{};

template<typename System, typename iterator>
struct parser_generator< typename TestModule1::type<System>::test_object1, System, iterator > {

    typedef rule<iterator, boost::shared_ptr<typename TestModule1::type<System>::test_object1>() > generator;
    static void init(generator& r) {
        r = string[karma::_1="<type>object 1 prop</type>\n<value>HaHAHAHAHA</value>"];
    };
};
}

typedef dcm::Kernel<double> Kernel;
typedef dcm::System<Kernel, dcm::ModuleState::type, TestModule1::type> System;

BOOST_AUTO_TEST_SUITE(parser_suit);
/*
BOOST_AUTO_TEST_CASE(parser_properties) {

    System sys;

    std::ostringstream s;
    dcm::ModuleState::type<System>::propperty_generator<test_prop> g(s);

    std::string test("value");
    g(test);

    std::cout<<s.str()<<std::endl;

    std::ostringstream s2;
    TestModule1::type<System>::test_object1 obj(sys);
    obj.setProperty<TestModule1::type<System>::test_object1_prop>(27);
    obj.setProperty<TestModule1::type<System>::test_object2_prop>("new value");

    dcm::ModuleState::type<System>::propperty_vector_generator<TestModule1::type<System>::test_object1::Sequence> functor(s2, obj.m_properties);
    mpl::for_each<TestModule1::type<System>::test_object1::Sequence>(functor);

    std::cout<<s2.str()<<std::endl;

}*/

BOOST_AUTO_TEST_CASE(parser_graph) {

    System sys;

    fusion::vector<dcm::LocalVertex, dcm::GlobalVertex> res1 = sys.m_cluster.addVertex();
    fusion::vector<dcm::LocalVertex, dcm::GlobalVertex> res2 = sys.m_cluster.addVertex();
    fusion::vector<dcm::LocalVertex, dcm::GlobalVertex> res3 = sys.m_cluster.addVertex();
    dcm::LocalEdge e1 = fusion::at_c<0>(sys.m_cluster.addEdge(fusion::at_c<1>(res1),fusion::at_c<1>(res2)));
    dcm::LocalEdge e2 = fusion::at_c<0>(sys.m_cluster.addEdge(fusion::at_c<1>(res2),fusion::at_c<1>(res3)));
    dcm::LocalEdge e3 = fusion::at_c<0>(sys.m_cluster.addEdge(fusion::at_c<1>(res1),fusion::at_c<1>(res3)));

    boost::shared_ptr<TestModule1::type<System>::test_object1> ptr(new TestModule1::type<System>::test_object1(sys));
    sys.m_cluster.setObject(fusion::at_c<1>(res1), ptr);
    
    sys.m_cluster.setProperty<TestModule1::type<System>::test_edge1_prop>(e1, 1);
    sys.m_cluster.setProperty<TestModule1::type<System>::test_edge1_prop>(e2, 2);
    sys.m_cluster.setProperty<TestModule1::type<System>::test_edge1_prop>(e3, 3);

    std::ostringstream s;
    sys.saveState(s);

    std::cout<<s.str()<<std::endl;
}

BOOST_AUTO_TEST_SUITE_END();
