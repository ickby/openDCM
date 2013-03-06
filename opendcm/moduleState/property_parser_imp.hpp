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

#ifndef DCM_PROPERTY_PARSER_IMP_H
#define DCM_PROPERTY_PARSER_IMP_H

#include "property_parser.hpp"


namespace dcm {

typedef boost::spirit::istream_iterator IIterator;

namespace details {

template<typename Prop>
empty_pop_parser<Prop>::empty_pop_parser(): empty_pop_parser<Prop>::base_type(start) {
    start = qi::eps(false);
};

template<typename Prop, typename Par>
prop_parser<Prop, Par>::prop_parser() : prop_parser<Prop, Par>::base_type(start) {
    Par::init(subrule);
    start =  qi::lit("<Property>") >> subrule >> qi::lit("</Property>");
};


template<typename Sys, typename PropertyList>
prop_par<Sys, PropertyList>::prop_par() : prop_par<Sys, PropertyList>::base_type(prop) {

    prop =       -(qi::eps(valid<0>::value) >> fusion::at<index<0> >(rules)[phx::at_c<index<0>::value>(qi::_val) = qi::_1])
                 >> -(qi::eps(valid<1>::value) >> fusion::at<index<1> >(rules)[phx::at_c<index<1>::value>(qi::_val) = qi::_1])
                 >> -(qi::eps(valid<2>::value) >> fusion::at<index<2> >(rules)[phx::at_c<index<2>::value>(qi::_val) = qi::_1])
                 >> -(qi::eps(valid<3>::value) >> fusion::at<index<3> >(rules)[phx::at_c<index<3>::value>(qi::_val) = qi::_1])
                 >> -(qi::eps(valid<4>::value) >> fusion::at<index<4> >(rules)[phx::at_c<index<4>::value>(qi::_val) = qi::_1])
                 >> -(qi::eps(valid<5>::value) >> fusion::at<index<5> >(rules)[phx::at_c<index<5>::value>(qi::_val) = qi::_1])
                 >> -(qi::eps(valid<6>::value) >> fusion::at<index<6> >(rules)[phx::at_c<index<6>::value>(qi::_val) = qi::_1])
                 >> -(qi::eps(valid<7>::value) >> fusion::at<index<7> >(rules)[phx::at_c<index<7>::value>(qi::_val) = qi::_1])
                 >> -(qi::eps(valid<8>::value) >> fusion::at<index<8> >(rules)[phx::at_c<index<8>::value>(qi::_val) = qi::_1])
                 >> -(qi::eps(valid<9>::value) >> fusion::at<index<9> >(rules)[phx::at_c<index<9>::value>(qi::_val) = qi::_1]);
};

template<typename Sys>
cluster_prop_par<Sys>::cluster_prop_par() : prop_par<Sys, typename Sys::Cluster::cluster_properties>() {};

template<typename Sys>
vertex_prop_par<Sys>::vertex_prop_par() : prop_par<Sys, typename Sys::Cluster::vertex_properties>() {};

template<typename Sys>
edge_prop_par<Sys>::edge_prop_par() : prop_par<Sys, typename Sys::Cluster::edge_properties>() {};

} //DCM
} //details

#endif
