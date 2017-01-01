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

#ifndef DCM_LOGGING_H
#define DCM_LOGGING_H

#ifdef DCM_USE_LOGGING

#define BOOST_LOG_DYN_LINK

#include <fstream>
#include <boost/log/core.hpp>
#include <boost/log/sinks.hpp>
#include <boost/log/expressions/formatters.hpp>
#include <boost/log/sources/severity_logger.hpp>
#include <boost/log/attributes.hpp>
#include <boost/log/sources/logger.hpp>
#include <boost/log/sources/record_ostream.hpp>
#include <boost/log/expressions.hpp>

#include <boost/shared_ptr.hpp>

namespace logging = boost::log;
namespace sinks = boost::log::sinks;
namespace src = boost::log::sources;
namespace expr = boost::log::expressions;
namespace attrs = boost::log::attributes;
namespace keywords = boost::log::keywords;

namespace dcm {
namespace details {

enum severity_level {

    iteration = 0,
    solving,
    manipulation,
    information,
    error
};

BOOST_LOG_ATTRIBUTE_KEYWORD(severity, "Severity", severity_level)

static int counter = 0;
typedef sinks::synchronous_sink< sinks::text_ostream_backend > sink_t;
typedef src::severity_logger< severity_level > dcm_logger;

inline boost::shared_ptr< sink_t > init_log() {

    //create the filename
    std::stringstream str;
    str<<"openDCM_"<<counter<<"_%N.log";
    counter++;

    typedef sinks::synchronous_sink< sinks::text_ostream_backend > text_sink;
    boost::shared_ptr< text_sink > sink = boost::make_shared< text_sink >();
    
    auto stream = boost::make_shared< std::ofstream >("openDCM.log");
    sink->locked_backend()->add_stream(stream);
   
    sink->set_formatter(
        expr::stream <<"[" << expr::attr<std::string>("Tag") <<"] "
        << expr::if_(expr::has_attr<std::string>("ID")) [
            expr::stream << "["<< expr::attr< std::string >("ID")<<"] "]
        << expr::smessage
    );

    logging::core::get()->add_sink(sink);
    return sink;
};

inline void stop_log(boost::shared_ptr< sink_t >& sink) {

    // Remove the sink from the core, so that no records are passed to it
    logging::core::get()->remove_sink(sink);
    sink.reset();
};

}; //details
}; //dcm

#endif //DCM_USE_LOGGING

#endif //DCM_LOGGING_H
