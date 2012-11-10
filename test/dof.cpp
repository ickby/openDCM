/*
    openDCM, dimensional constraint manager
    Copyright (C) 2012  Stefan Troeger <stefantroeger@gmx.net>

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

#include "opendcm/core/kernel.hpp"
#include "opendcm/module3d/dof.hpp"

#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(dof);

BOOST_AUTO_TEST_CASE(dof_translation) {

    typedef dcm::Kernel<double> Kernel;
    typedef typename dcm::Dof<Kernel, int>::ConstraintVector::iterator iter;
    typedef dcm::Dof<Kernel, int>::Result Result;
    dcm::Dof<Kernel, int> d;

    BOOST_CHECK(d.dofTranslational() == 3);
    BOOST_CHECK(d.dofRotational() == 3);
    BOOST_CHECK(d.dof() == 6);

    Kernel::Vector3 v(1,0,0);
    Result r = d.removeTranslationDirection(v, 1);
    BOOST_CHECK(r.first);

    BOOST_CHECK(d.dofTranslational() == 2);
    BOOST_CHECK(d.dofRotational() == 3);
    BOOST_CHECK(d.dof() == 5);

    r = d.removeTranslationDirection(v,2);
    BOOST_CHECK(d.dofTranslational() == 2);
    BOOST_CHECK(!r.first);
    BOOST_CHECK(r.second.front() == 1);
    BOOST_CHECK(r.second.size() == 1);

    Kernel::Vector3 v2(0,1,0);
    r = d.removeTranslationDirection(v2,3);
    BOOST_CHECK(d.dofTranslational() == 1);
    BOOST_CHECK(r.first);

    Kernel::Vector3 v3(1,1,0);
    r = d.removeTranslationDirection(v3,4);
    BOOST_CHECK(d.dofTranslational() == 1);
    BOOST_CHECK(!r.first);
    iter it = r.second.begin();
    BOOST_CHECK(*it == 1);
    BOOST_CHECK(*(++it) == 3);
    BOOST_CHECK(++it == r.second.end());

    Kernel::Vector3 v4(1,1,1);
    r = d.removeTranslationDirection(v4,5);
    BOOST_CHECK(d.dofTranslational() == 0);
    BOOST_CHECK(r.first);

    Kernel::Vector3 v5(7,2,5);
    r = d.removeTranslationDirection(v5,6);
    BOOST_CHECK(d.dofTranslational() == 0);
    BOOST_CHECK(!r.first);
    it = r.second.begin();
    BOOST_CHECK(*it == 1);
    BOOST_CHECK(*(++it) == 3);
    BOOST_CHECK(*(++it) == 5);
    BOOST_CHECK(++it == r.second.end());

};

BOOST_AUTO_TEST_CASE(dof_rotational) {

    typedef dcm::Kernel<double> Kernel;
    typedef typename dcm::Dof<Kernel, int>::ConstraintVector::iterator iter;
    typedef dcm::Dof<Kernel, int>::Result Result;
    dcm::Dof<Kernel, int> d;

    Kernel::Vector3 v(1,0,0);
    Result r = d.allowOnlyRotationDirection(v, 1);
    BOOST_CHECK(r.first);

    BOOST_CHECK(d.dofTranslational() == 3);
    BOOST_CHECK(d.dofRotational() == 1);
    BOOST_CHECK(d.dof() == 4);

    r = d.allowOnlyRotationDirection(v,2);
    BOOST_CHECK(d.dofRotational() == 1);
    BOOST_CHECK(!r.first);
    BOOST_CHECK(r.second.front() == 1);
    BOOST_CHECK(r.second.size() == 1);

    Kernel::Vector3 v2(0,1,0);
    r = d.allowOnlyRotationDirection(v2,3);
    BOOST_CHECK(d.dofRotational() == 0);
    BOOST_CHECK(r.first);

    Kernel::Vector3 v5(7,2,5);
    r = d.allowOnlyRotationDirection(v5,6);
    BOOST_CHECK(d.dofRotational() == 0);
    BOOST_CHECK(!r.first);
    iter it = r.second.begin();
    BOOST_CHECK(*it == 1);
    BOOST_CHECK(*(++it) == 3);
    BOOST_CHECK(++it == r.second.end());

};

BOOST_AUTO_TEST_SUITE_END();
