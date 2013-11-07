/*
    openDCM, dimensional constraint manager
    Copyright (C) 2012  Stefan Troeger <stefantroeger@gmx.net>

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

#ifndef DCM_MODULE3D_H
#define DCM_MODULE3D_H

#define DCM_USE_MODULE3D

#ifdef _WIN32
	//warning about to long decoraded names, won't affect the code correctness
	#pragma warning( disable : 4503 )
#endif

#include "module3d/defines.hpp"
#include "module3d/geometry.hpp"
#include "module3d/distance.hpp"
#include "module3d/parallel.hpp"
#include "module3d/angle.hpp"
#include "module3d/coincident.hpp"
#include "module3d/alignment.hpp"
#include "module3d/module.hpp"

#ifdef DCM_USE_MODULESTATE
#include "module3d/state.hpp"

#ifdef USE_EXTERNAL
#define DCM_EXTERNAL_INCLUDE_009 <opendcm/module3d/state_imp.hpp>
#define DCM_EXTERNAL_009( System )\
    template struct dcm::parser_generator< typename dcm::details::getModule3D< System >::type::Geometry3D , System, std::ostream_iterator<char> >; \
    template struct dcm::parser_generator< typename dcm::details::getModule3D<System>::type::vertex_prop , System, std::ostream_iterator<char> >; \
    template struct dcm::parser_generator< typename dcm::details::getModule3D<System>::type::Constraint3D , System, std::ostream_iterator<char> >; \
    template struct dcm::parser_generator< typename dcm::details::getModule3D<System>::type::edge_prop , System, std::ostream_iterator<char> >; \
    template struct dcm::parser_generator< typename dcm::details::getModule3D<System>::type::fix_prop, System, std::ostream_iterator<char> >; \
    template struct dcm::parser_parser< typename dcm::details::getModule3D<System>::type::edge_prop, System, boost::spirit::istream_iterator >; \
    template struct dcm::parser_parser< typename dcm::details::getModule3D<System>::type::vertex_prop, System, boost::spirit::istream_iterator >; \
    template struct dcm::parser_parser< typename dcm::details::getModule3D<System>::type::fix_prop, System, boost::spirit::istream_iterator >; \
    template struct dcm::parser_parser< typename dcm::details::getModule3D<System>::type::Geometry3D, System, boost::spirit::istream_iterator >;

#define DCM_EXTERNAL_INCLUDE_010 <opendcm/module3d/state_imp.hpp>
#define DCM_EXTERNAL_010( System )\
    template struct dcm::parser_parser< typename dcm::details::getModule3D<System>::type::Constraint3D, System, boost::spirit::istream_iterator >;

    
#endif //ecternal
#endif //use state

#endif //DCM_MODULE3D_H

