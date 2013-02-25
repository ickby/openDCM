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

#ifndef DCM_EXTRACTOR_H
#define DCM_EXTRACTOR_H

#include <opendcm/core/clustergraph.hpp>
#include <boost/fusion/include/at_c.hpp>

namespace fusion = boost::fusion;

namespace dcm {
  
typedef std::ostream_iterator<char> Iterator;

template<typename Sys>
struct Extractor {
  
    typedef typename boost::graph_traits<typename Sys::Cluster>::vertex_iterator viter;
    typedef typename boost::graph_traits<typename Sys::Cluster>::edge_iterator eiter;
  
    void getVertexRange(typename Sys::Cluster* cluster, std::vector<typename Sys::Cluster::vertex_bundle>& range) {
        std::pair<viter, viter> res = boost::vertices(*cluster);
        for(; res.first != res.second; res.first++)
            range.push_back((*cluster)[*res.first]);
    };
    void getEdgeRange(typename Sys::Cluster* cluster,
                      std::vector<fusion::vector3<typename Sys::Cluster::edge_bundle, GlobalVertex, GlobalVertex> >& range) {

        std::pair<eiter, eiter> res = boost::edges(*cluster);
        for(; res.first != res.second; res.first++)
            range.push_back(fusion::make_vector((*cluster)[*res.first],
                                                (*cluster).getGlobalVertex(boost::source(*res.first, *cluster)),
                                                (*cluster).getGlobalVertex(boost::target(*res.first, *cluster))));

    };
    void getGlobalEdgeSource(typename Sys::Cluster::edge_bundle_single b, int& source) {
        source = fusion::at_c<1>(b).source;
    };
    void getGlobalEdgeTarget(typename Sys::Cluster::edge_bundle_single b, int& target) {
        target = fusion::at_c<1>(b).target;
    };
    void getGlobalEdgeID(typename Sys::Cluster::edge_bundle_single b, int& id) {
        id = fusion::at_c<1>(b).ID;
    };
    void getVertexID(typename Sys::Cluster* cluster, LocalVertex v, int& id) {
        if(v)
            id = cluster->getGlobalVertex(v);
        else id=0;
    };
    void makeInitPair(typename Sys::Cluster& cluster, std::pair<LocalVertex, typename Sys::Cluster*>& pair) {
        pair = std::make_pair((LocalVertex)NULL, &cluster);
    };
    void getClusterRange(typename Sys::Cluster* cluster, std::map<LocalVertex, typename Sys::Cluster*>& range) {
        range = cluster->m_clusters;
    };
    void getClusterProperties(typename Sys::Cluster* cluster, typename Sys::Cluster::cluster_bundle& seq) {
        seq = cluster->m_cluster_bundle;
    };
};



}//namespace dcm

#endif //DCM_GENERATOR_H




