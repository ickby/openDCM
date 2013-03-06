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

#ifndef DCM_PARSER_H
#define DCM_PARSER_H

#include <iosfwd>

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/support_istream_iterator.hpp>
#include <boost/spirit/include/qi_string.hpp>
#include <boost/spirit/include/phoenix.hpp>

#include "opendcm/core/clustergraph.hpp"

#include "property_parser.hpp"
#include "object_parser.hpp"
#include "edge_vertex_parser.hpp"
#include "extractor.hpp"

namespace qi = boost::spirit::qi;
namespace ascii = boost::spirit::ascii;
namespace phx = boost::phoenix;

namespace dcm {

typedef boost::spirit::istream_iterator IIterator;

struct sp : qi::grammar<IIterator, std::string()> {

    qi::rule<IIterator, std::string()> start;
    sp() : sp::base_type(start) {
    start %= +qi::char_;
};
static void print(std::string& s) {
    std::cout<<"parsed string:"<<std::endl<<s<<std::endl<<"done print string"<<std::endl;
};
};

template<typename Sys>
struct parser : qi::grammar<IIterator, Sys*(), qi::space_type> {

    parser(Sys& s);

    qi::rule<IIterator, Sys*(), qi::space_type> start;

    qi::rule<IIterator, fusion::vector2<GlobalVertex, typename Sys::Cluster*>(Sys*), qi::space_type> cluster;
    details::cluster_prop_par<Sys> cluster_prop;
    
    details::obj_par<Sys> objects;
    
    details::vertex_parser<Sys> vertex;
    details::edge_parser<Sys> edge;

    sp str;
    Sys& system;
    Injector<Sys> in;
};

}

#ifndef USE_EXTERNAL
#include "parser_imp.hpp"
#endif

#endif //DCM_PARSER_H
