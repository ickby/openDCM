/*
    openDCM, dimensional constraint manager
    Copyright (C) 2012  Stefan Troeger <stefantroeger@gmx.net>

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

#ifndef GCM_LOGGING_H
#define GCM_LOGGING_H

#ifdef USE_LOGGING

#include <boost/log/core.hpp>
#include <boost/log/sinks.hpp>
#include <boost/log/attributes.hpp>
#include <boost/log/attributes/attribute_def.hpp>
#include <boost/log/utility/init/to_file.hpp>
#include <boost/log/utility/init/common_attributes.hpp>
#include <boost/log/sources/logger.hpp>
#include <boost/log/sources/record_ostream.hpp>
#include <boost/log/formatters.hpp>
#include <boost/log/filters.hpp>


namespace logging = boost::log;
//namespace sinks = boost::log::sinks;
namespace src = boost::log::sources;
namespace fmt = boost::log::formatters;
namespace flt = boost::log::filters;
namespace attrs = boost::log::attributes;
namespace keywords = boost::log::keywords;

namespace dcm {
void init_log() {
  
    //regiser a sink for file ouput
    logging::init_log_to_file(
	keywords::file_name = "openDCM_%N.log",                 // file name pattern
        keywords::rotation_size = 10 * 1024 * 1024, 		//Rotation size is 10MB
	keywords::format =
        (
           fmt::stream <<"[" << fmt::attr<std::string>("Tag") <<"] " << fmt::message()
        )
    );
    //logging::add_common_attributes();
};
} //dcm

#endif //USE_LOGGING

#endif //GCM_LOGGING_H
