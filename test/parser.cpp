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
#include "opendcm/ModuleParser"

#include <boost/test/unit_test.hpp>

struct TestModule1 {

    template<typename Sys>
    struct type {
        typedef mpl::map0<> signal_map;

        struct test_object1 : public dcm::Object<Sys, test_object1, signal_map > {
            test_object1(Sys& system) : dcm::Object<Sys, test_object1, signal_map >(system) { };
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

        typedef mpl::vector1<test_object1> objects;
        typedef mpl::vector3<test_object1_prop>   properties;

        static void system_init(Sys& sys) {};
    };
};

BOOST_AUTO_TEST_SUITE(parser_suit);

BOOST_AUTO_TEST_CASE(parser_properties) {

    

}

BOOST_AUTO_TEST_SUITE_END();