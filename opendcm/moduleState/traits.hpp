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

#ifndef DCM_PARSER_TRAITS_H
#define DCM_PARSER_TRAITS_H

#include <boost/mpl/bool.hpp>

#include <boost/spirit/include/karma.hpp>
#include <boost/spirit/include/karma_string.hpp>
#include <boost/spirit/include/karma_rule.hpp>

namespace dcm {
  
using boost::spirit::karma::rule;
using boost::spirit::ascii::string;

template<typename type, typename System>
struct parser_generate : public mpl::false_ {};

template<typename type, typename System, typename iterator>
struct parser_generator {
  typedef rule<iterator> generator;
  
  static void init(generator& r) {
    r = boost::spirit::karma::eps(false);
  };
};
}
#endif //DCM_PARSER_TRAITS_H
