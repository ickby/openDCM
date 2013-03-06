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

#ifndef DCM_PARSER_IMP_H
#define DCM_PARSER_IMP_H

#include "parser.hpp"

namespace dcm {

typedef boost::spirit::istream_iterator IIterator;

template<typename Sys>
parser<Sys>::parser(Sys& s) : parser<Sys>::base_type(start), system(s) {


    cluster = qi::lit("<Cluster id=")[phx::at_c<1>(qi::_val) = phx::new_<typename Sys::Cluster>()]
              >> qi::int_[phx::at_c<0>(qi::_val) = qi::_1] >> ">"
              >> cluster_prop[phx::bind(&Injector<Sys>::setClusterProperties, in, phx::at_c<1>(qi::_val), qi::_1)]
              >> *vertex(phx::at_c<1>(qi::_val), qi::_r1)
              >> *edge(phx::at_c<1>(qi::_val), qi::_r1)
              >> *(cluster(qi::_r1)[phx::bind(&Injector<Sys>::addCluster, in, phx::at_c<1>(qi::_val), qi::_1)])
              >> "</Cluster>";

    start = qi::no_skip[ascii::string("<Cluster id=0>")] >> cluster_prop[phx::bind(&Injector<Sys>::setClusterProperties, in, &phx::bind(&Sys::m_cluster, qi::_val), qi::_1)]
            >> *(vertex(&phx::bind(&Sys::m_cluster, qi::_val), qi::_val))
            >> *(edge(&phx::bind(&Sys::m_cluster, qi::_val), qi::_val))
            >> *(cluster(qi::_val)[phx::bind(&Injector<Sys>::addCluster, in, &phx::bind(&Sys::m_cluster, qi::_val), qi::_1)])
            >> "</Cluster>" >> str[phx::bind(&sp::print, qi::_1)];
};

}
#endif //DCM_PARSER_H
