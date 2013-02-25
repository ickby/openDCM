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

#ifndef DCM_EDGE_GENERATOR_IMP_H
#define DCM_EDGE_GENERATOR_IMP_H

#include "edge_vertex_generator.hpp"

namespace dcm {
namespace details {
  
template<typename Sys>
edge_generator<Sys>::edge_generator() : edge_generator<Sys>::base_type(edge_range) {
  
        globaledge = int_[phx::bind(&Extractor<Sys>::getGlobalEdgeID, ex, _val, karma::_1)]
                     << " source=" << int_[phx::bind(&Extractor<Sys>::getGlobalEdgeSource, ex, _val, karma::_1)]
                     << " target=" << int_[phx::bind(&Extractor<Sys>::getGlobalEdgeTarget, ex, _val, karma::_1)] << '>'
                     << "+" << objects[karma::_1 = phx::at_c<0>(_val)] << "-\n" ;


        globaledge_range = *(lit("<GlobalEdge id=")<<globaledge<<lit("</GlobalEdge>"));

        edge =  lit("source=")<<int_[karma::_1 = phx::at_c<1>(_val)] << " target="<<int_[karma::_1 = phx::at_c<2>(_val)] << ">+"
                << edge_prop[karma::_1 = phx::at_c<0>(phx::at_c<0>(_val))]
                << eol << globaledge_range[karma::_1 = phx::at_c<1>(phx::at_c<0>(_val))] << '-' << eol;

        edge_range = (lit("<Edge ") << edge << lit("</Edge>")) % eol;
};

template<typename Sys>
vertex_generator<Sys>::vertex_generator() : vertex_generator<Sys>::base_type(vertex_range) {
  
        vertex = int_[karma::_1 = phx::at_c<2>(_val)] << ">+"
                 << vertex_prop[karma::_1 = phx::at_c<0>(_val)]
                 << objects[karma::_1 = phx::at_c<1>(_val)]
                 << "-\n";

        vertex_range = '\n' << (lit("<Vertex id=") << vertex  << lit("</Vertex>")) % eol;
};

}//details
}//dcm

#endif