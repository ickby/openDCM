/*
    openDCM, dimensional constraint manager
    Copyright (C) 2014  Stefan Troeger <stefantroeger@gmx.net>

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

#ifndef DCM_MODULE_CORE_H
#define DCM_MODULE_CORE_H

#include <boost/mpl/int.hpp>

#include "clustergraph.hpp"
#include "logging.hpp"
#include "object.hpp"

namespace dcm {
namespace details {

namespace mpl = boost::mpl;

template<typename Final>
struct ModuleCoreInit {

    ModuleCoreInit() : graph(NULL)
#ifdef USE_LOGGING
        , sink(init_log())
#endif
    {};

    ~ModuleCoreInit() {
        if(graph)
            delete graph;

#ifdef USE_LOGGING
        stop_log(sink);
#endif
    };

#ifdef USE_LOGGING
    template<typename Expr>
    void setLoggingFilter(const Expr& ex) {
        sink->set_filter(ex);
    }
#endif

    //initialize the objectbase for unified object handling
    typedef Object<Final> ObjectBase;

protected:
    typedef mpl::vector0<> EdgeProperties;
    typedef mpl::vector0<> GlobalEdgeProperties;
    typedef mpl::vector0<> VertexProperties;
    typedef mpl::vector0<> ClusterProperties;

    //ensure that the correct graph type is used by not allowing anyone to set the graph pointer
    ClusterGraphBase* getGraph() {
        if(!graph)
            graph = new typename Final::Graph();

        return graph;
    };

private:
    ClusterGraphBase* graph;
#ifdef USE_LOGGING
    boost::shared_ptr< sink_t > sink;
#endif

};

template<typename Final, typename Stacked>
struct ModuleCoreFinish : public Stacked {

protected:
    typedef ClusterGraph<typename Stacked::EdgeProperties, typename Stacked::VertexProperties,
            typename Stacked::ClusterProperties, typename Stacked::GlobalEdgeProperties> Graph;
};

} //details
} //dcm
#endif //DCM_MODULE_CORE_H
