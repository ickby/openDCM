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

//move the traits specializations outside of the traits definition to avoid the spirit header parsing every
//time this module is included and just parse it in externalisation mode when the generator is build

#ifndef DCM_PARSER_TRAITS_IMPL_H
#define DCM_PARSER_TRAITS_IMPL_H

#include "traits.hpp"

#include <boost/mpl/bool.hpp>

#include <boost/spirit/include/karma.hpp>
#include <boost/spirit/include/karma_string.hpp>
#include <boost/spirit/include/karma_int.hpp>
#include <boost/spirit/include/karma_bool.hpp>
#include <boost/spirit/include/karma_rule.hpp>
#include <boost/spirit/include/karma_auto.hpp>

namespace karma = boost::spirit::karma;

namespace boost {
namespace spirit {
namespace traits {
template <>
struct create_generator<dcm::No_Identifier> {

    typedef BOOST_TYPEOF(karma::eps(false)) type;
    static type call()  {
        return karma::eps(false);
    }
};
}
}
}

namespace dcm {
  
  template<typename System>
struct parser_generate<type_prop, System> : public mpl::true_ {};

template<typename System, typename iterator>
struct parser_generator<type_prop, System, iterator> {
    typedef karma::rule<iterator, int()> generator;

    static void init(generator& r) {
        r = karma::lit("<type>clustertype</type>\n<value>") << karma::int_ <<"</value>";
    };
};

template<typename System>
struct parser_generate<changed_prop, System> : public mpl::true_ {};

template<typename System, typename iterator>
struct parser_generator<changed_prop, System, iterator> {
    typedef karma::rule<iterator, bool()> generator;

    static void init(generator& r) {
        r = karma::lit("<type>clusterchanged</type>\n<value>") << karma::bool_ <<"</value>";
    };
};

template<typename System>
struct parser_generate<id_prop<typename System::Identifier>, System>
        : public mpl::not_<boost::is_same<typename System::Identifier, No_Identifier> > {};

template<typename System, typename iterator>
struct parser_generator<id_prop<typename System::Identifier>, System, iterator> {
    typedef karma::rule<iterator, typename System::Identifier()> generator;

    static void init(generator& r) {
        r = karma::lit("<type>id</type>\n<value>") << karma::auto_ <<"</value>";
    };
};

} //namespace dcm

#endif //DCM_PARSER_TRAITS_IMPL_H