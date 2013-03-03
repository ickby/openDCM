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

#ifndef DCM_EXTERNALIZE_H
#define DCM_EXTERNALIZE_H

#include <boost/preprocessor/list/for_each.hpp>
#include <boost/preprocessor/cat.hpp>
#include <boost/preprocessor/list/append.hpp>

#ifdef USE_EXTERNAL

#define FOLD(r, data, elem) extern template class elem<data>;

#ifdef DCM_USE_MODULESTATE

#define LIST_STATE (dcm::generator, \
                    (dcm::details::obj_gen, \
                     (dcm::details::vertex_prop_gen, \
                      (dcm::details::edge_prop_gen, \
                       (dcm::details::cluster_prop_gen, \
                        (dcm::details::edge_generator, \
                         (dcm::details::vertex_generator, BOOST_PP_NIL)))))))                     
#else
#define LIST_STATE BOOST_PP_NIL
#endif

#ifdef DCM_USE_MODULE3D
#define LIST_3D BOOST_PP_NIL
#else
#define LIST_3D BOOST_PP_NIL
#endif

#define LIST BOOST_PP_LIST_APPEND(LIST_STATE, LIST_3D)

#define DCM_EXTERNALIZE( System ) \
    BOOST_PP_LIST_FOR_EACH(FOLD,System, LIST)

#define DCM_EXTERNAL_INCLUDE_001 <opendcm/moduleState/generator.hpp>
#define DCM_EXTERNAL_001( System )\
    template class dcm::generator<System>; \
     
#define DCM_EXTERNAL_INCLUDE_002 <opendcm/moduleState/object_generator_imp.hpp>
#define DCM_EXTERNAL_002( System )\
    template class dcm::details::obj_gen<System>; \
     
#define DCM_EXTERNAL_INCLUDE_003 <opendcm/moduleState/property_generator_imp.hpp>
#define DCM_EXTERNAL_003( System )\
    template class dcm::details::vertex_prop_gen<System>; \
    template class dcm::details::edge_prop_gen<System>; \
    template class dcm::details::cluster_prop_gen<System>;

#define DCM_EXTERNAL_INCLUDE_004 <opendcm/moduleState/edge_vertex_generator_imp.hpp>
#define DCM_EXTERNAL_004( System )\
    template class dcm::details::edge_generator<System>; \
    template class dcm::details::vertex_generator<System>; \
     
#else
#define DCM_EXTERNALIZE( System )
#endif

#endif //DCM_EXTERNALIZE_H

