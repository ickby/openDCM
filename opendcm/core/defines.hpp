/*
    openDCM, dimensional constraint manager
    Copyright (C) 2013  Stefan Troeger <stefantroeger@gmx.net>

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

#ifndef DCM_DEFINES_CORE_H
#define DCM_DEFINES_CORE_H


#include <boost/exception/exception.hpp>
#include <boost/graph/graph_concepts.hpp>
#include <string>

namespace dcm {

//all solving related errors
typedef boost::error_info<struct user_message,std::string> error_message;
typedef boost::error_info<struct first_geom, std::string>  error_type_first_geometry;
typedef boost::error_info<struct second_geom, std::string> error_type_second_geometry;

//exception codes are needed by the user
struct creation_error : virtual boost::exception {};
struct solving_error : virtual boost::exception {};
struct constraint_error : virtual boost::exception {};
struct general_assert_error : virtual boost::exception {};

//we may need to throw exceptions instead of failing the application in assert as we only can detect
//exceptions in the testing framework. Therefore we write our own assert
#ifdef DCM_DEBUG 
#ifdef DCM_DEBUG_THROW_AT_ASSERT
#define dcm_assert(exp) \
if(exp) \
    throw general_assert_error();
#else
#define dcm_assert(exp) \
    assert(exp);
#endif
#else
#define dcm_assert(exp)
#endif
}; //dcm

#endif //DCM_DEFINES_CORE_H
