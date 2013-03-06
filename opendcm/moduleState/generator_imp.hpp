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

#ifndef DCM_GENERATOR_IMP_H
#define DCM_GENERATOR_IMP_H

#include "generator.hpp"

namespace dcm {

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



