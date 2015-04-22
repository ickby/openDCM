/*
    openDCM, dimensional constraint manager
    Copyright (C) 2015  Stefan Troeger <stefantroeger@gmx.net>

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

#include "opendcm/core/property.hpp"

#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(property);

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

TestContainer1 : public dcm::details::PropertyContainer {
    
    dcm::details::Property<TestProperty1> testProperty1;
};

BOOST_AUTO_TEST_CASE(basics) {

    System sys;

    BOOST_CHECK(sys.module_function1() == 1);
    BOOST_CHECK(sys.module_function2() == 2);

};

BOOST_AUTO_TEST_SUITE_END();

