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

#ifndef DCM_GENERATOR_H
#define DCM_GENERATOR_H

#include "object_generator.hpp"
#include "property_generator.hpp"
#include "extractor.hpp"

#include <opendcm/core/clustergraph.hpp>

#include "traits.hpp"
#include "traits_impl.hpp"
#include "indent.hpp"
#include "edge_vertex_generator.hpp"

#include <boost/mpl/int.hpp>
#include <boost/mpl/for_each.hpp>
#include <boost/fusion/container/vector/convert.hpp>
#include <boost/fusion/include/as_vector.hpp>
#include <boost/fusion/include/at.hpp>
#include <boost/fusion/include/at_c.hpp>
#include <boost/fusion/include/std_pair.hpp>

#include <boost/spirit/include/karma.hpp>
#include <boost/spirit/include/phoenix.hpp>

namespace karma = boost::spirit::karma;
namespace phx = boost::phoenix;
namespace fusion = boost::fusion;

namespace dcm {

typedef std::ostream_iterator<char> Iterator;

template<typename Sys>
struct generator : karma::grammar<Iterator, typename Sys::Cluster& ()> {

    typedef typename boost::graph_traits<typename Sys::Cluster>::vertex_iterator viter;
    typedef typename boost::graph_traits<typename Sys::Cluster>::edge_iterator eiter;

    generator(Sys& s);

    karma::rule<Iterator, typename Sys::Cluster& ()> start;

    karma::rule<Iterator, std::map<LocalVertex, typename Sys::Cluster*>()> cluster_range;
    karma::rule<Iterator, std::pair<LocalVertex, typename Sys::Cluster*>()> cluster;
    details::cluster_prop_gen<Sys> cluster_prop;

    details::vertex_generator<Sys> vertex_range;
    details::edge_generator<Sys> edge_range;

    details::obj_gen<Sys> objects;
    Extractor<Sys> ex;

    Sys& system;
};

//definitions outside to avoid automatic inline
template<typename Sys>
generator<Sys>::generator(Sys& s) : generator<Sys>::base_type(start), system(s) {

    cluster = karma::lit("<Cluster id=") <<karma::int_[phx::bind(&Extractor<Sys>::getVertexID, ex, phx::at_c<1>(karma::_val), phx::at_c<0>(karma::_val), karma::_1)]
              << ">+" << cluster_prop[phx::bind(&Extractor<Sys>::getClusterProperties, ex, phx::at_c<1>(karma::_val), karma::_1)]
              << vertex_range[phx::bind(&Extractor<Sys>::getVertexRange, ex, phx::at_c<1>(karma::_val), karma::_1)]
              << -karma::buffer["\n" << edge_range[phx::bind(&Extractor<Sys>::getEdgeRange, ex, phx::at_c<1>(karma::_val), karma::_1)]]
              << -karma::buffer["\n" << cluster_range[phx::bind(&Extractor<Sys>::getClusterRange, ex, phx::at_c<1>(karma::_val), karma::_1)]] << "-\n"
              << "</Cluster>";

    cluster_range = cluster % karma::eol;

    start = cluster[phx::bind(&Extractor<Sys>::makeInitPair, ex, karma::_val, karma::_1)];
};

}//namespace dcm

#endif //DCM_GENERATOR_H



