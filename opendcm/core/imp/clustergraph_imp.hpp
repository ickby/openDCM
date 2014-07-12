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

#ifndef DCM_CLUSTERGRAPH_IMP_HPP
#define DCM_CLUSTERGRAPH_IMP_HPP

#include "../clustergraph.hpp"

#include <boost/fusion/include/algorithm.hpp>
#include <boost/fusion/include/at.hpp>
#include <boost/fusion/include/at_c.hpp>
#include <boost/fusion/include/make_vector.hpp>

#include <boost/bind.hpp>

namespace dcm {


//functors needed for implementation only
//***************************************

template<typename Graph>
struct edge_copier {
    edge_copier(const Graph& g1, Graph& g2)
        : graph1(g1), graph2(g2) { }

    void operator()(LocalEdge e1, LocalEdge e2) const {
        graph2[e2] = graph1[e1];
    }
    const Graph& graph1;
    Graph& graph2;
};

template<typename Graph>
struct vertex_copier {
    vertex_copier(const Graph& g1, Graph& g2)
        : graph1(g1), graph2(g2) { }

    void operator()(LocalVertex v1, LocalVertex v2) const {
        graph2[v2] = graph1[v1];
    }
    const Graph& graph1;
    Graph& graph2;
};

struct placehoder {
    template<typename T>
    void operator()(T t) {};
};

template<typename Functor, typename Graph>
struct apply_remove_prediacte {
    Functor& func;
    GlobalVertex vert;
    GlobalEdge edge;
    bool isEdge;

    apply_remove_prediacte(Functor& f, GlobalVertex v) : func(f), vert(v), isEdge(false) {};
    apply_remove_prediacte(Functor& f, GlobalEdge e) : func(f), edge(e), vert(0), isEdge(true) {};
    bool operator()(typename Graph::edge_bundle_single& e) {
        bool res;

        //this predicate can be used to compare the edge itself or the vertives it connects. See
        //if we are a relevant edge
        if(isEdge)
            res = (edge == fusion::at_c<0> (e));
        else
            res = (vert == fusion::at_c<0> (e).source) || (vert == fusion::at_c<0> (e).target);

        //we are a hit, invoke the functor.
        if(res || vert < 0)
            func(fusion::at_c<0> (e));

        return res || vert < 0;
    }
};

template<typename prop, typename Graph>
struct get_property_helper {

    typedef typename prop::type base_type;
    typedef const base_type& result_type;
    typedef typename mpl::if_< is_edge_property<prop, Graph>, mpl::int_<0>, mpl::int_<1> >::type pos;

    template<typename bundle>
    result_type operator()(bundle& p) {
        return fusion::at<pos>(p).template getProperty<prop>();
    }
};

template<typename prop, typename Graph>
struct set_property_helper {

    typedef typename prop::type base_type;
    typedef const base_type& value_type;
    typedef void result_type;
    typedef typename mpl::if_< is_edge_property<prop, Graph>, mpl::int_<0>, mpl::int_<1> >::type pos;

    set_property_helper(value_type v) :  value(v) {};

    template<typename bundle>
    result_type operator()(bundle& p) {
        fusion::at<pos>(p).template setProperty<prop>(value);
    };

    value_type value;
};

template<typename prop, typename Graph>
struct change_property_helper {

    typedef bool result_type;
    typedef typename mpl::if_< is_edge_property<prop, Graph>, mpl::int_<0>, mpl::int_<1> >::type pos;

    template<typename bundle>
    result_type operator()(bundle& p) {
        return fusion::at<pos>(p).template isPropertyChanged<prop>();
    }
};

template<typename prop, typename Graph>
struct acknowledge_property_helper {

    typedef void result_type;
    typedef typename mpl::if_< is_edge_property<prop, Graph>, mpl::int_<0>, mpl::int_<1> >::type pos;

    template<typename bundle>
    result_type operator()(bundle& p) {
        fusion::at<pos>(p).template acknowledgePropertyChange<prop>();
    }
};

template<typename Graph, typename Key>
struct changes_property_helper {

    typedef bool result_type;
    typedef typename mpl::if_< boost::is_same<Key, LocalEdge>, mpl::int_<0>, mpl::int_<1> >::type pos;

    template<typename bundle>
    result_type operator()(bundle& p) {
        return fusion::at<pos>(p).template hasPropertyChanges();
    }
};

template<typename Graph, typename Key>
struct ack_changes_property_helper {

    typedef void result_type;
    typedef typename mpl::if_< boost::is_same<Key, LocalEdge>, mpl::int_<0>, mpl::int_<1> >::type pos;

    template<typename bundle>
    result_type operator()(bundle& p) {
        fusion::at<pos>(p).template acknowledgePropertyChanges();
    }
};

//Function implementation
//***********************

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
template<typename T>
typename ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::global_extractor::result_type
ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::global_extractor::operator()(T& bundle) const {
    return fusion::at_c<0> (bundle);
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::global_vertex_extractor::global_vertex_extractor(ClusterGraph& g) : graph(g) {};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
typename ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::global_vertex_extractor::result_type
ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::global_vertex_extractor::operator()(LocalVertex& v) const {
    return graph.getGlobalVertex(v);
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
template<typename Functor>
void ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::copyInto(boost::shared_ptr<ClusterGraph> into, Functor& functor) const {

    //lists does not provide vertex index, so we have to build our own (cant use the internal
    //vertex_indexerty as we would need to reset the indices and that's not possible in const graph)
    typedef std::map<LocalVertex, int> IndexMap;
    IndexMap mapIndex;
    boost::associative_property_map<IndexMap> propmapIndex(mapIndex);

    std::pair<local_vertex_iterator, local_vertex_iterator>  vit = boost::vertices(*this);

    for(int c = 0; vit.first != vit.second; vit.first++, c++)
        put(propmapIndex, *vit.first, c);

    //first copy all vertices and edges, but be aware that the objects in the new graph
    //are also copys only and point to the old graph. there is a bug in older boost version
    //(<1.5 i belive) that breaks vertex_all propety map for bundled properties, so we
    //have to create our own copie functors
    into->clear();
    vertex_copier<Graph> vc(*this, *into);
    edge_copier<Graph> ec(*this, *into);
    boost::copy_graph(*this, *into, boost::vertex_index_map(propmapIndex).vertex_copy(vc).edge_copy(ec));

    //set the IDgen to the same value to avoid duplicate id's in the copied cluster
    into->m_id->setCount(m_id->count());

    //now that we have all vertices we can recreate the subclusters
    std::pair<const_cluster_iterator, const_cluster_iterator> it = clusters();

    for(; it.first != it.second; it.first++) {
        //create the new Graph
        boost::shared_ptr<ClusterGraph> ng = boost::shared_ptr<ClusterGraph> (new ClusterGraph(into));

        //we already have the new vertex, however, we need to find it
        GlobalVertex gv = getGlobalVertex((*it.first).first);
        LocalVertex  lv = into->getLocalVertex(gv).first;

        //add the new graph to the subclustermap
        into->m_clusters[lv] = ng;

        //copy the subcluster
        (*it.first).second->copyInto(ng, functor);
    }

    //lets see if the objects need special treatment
    into->for_each_object(functor, false);
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
template<typename T>
bool ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::operator== (const T& other) const {
    return this == &other;
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
template<typename T>
bool ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::operator!= (const T& other) const {
    return !(this == &other);
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
void ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::setCopyMode(bool on) {
    copy_mode = on;
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
template<typename P>
typename P::type& ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::getSubclusterProperty(LocalVertex v) {
    return getVertexCluster(v)->template getProperty<P>();
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
std::pair<boost::shared_ptr< ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop> >, LocalVertex> ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::createCluster() {
    vertex_bundle vp;
    fusion::at_c<0> (vp) = m_id->generate();
    LocalVertex v = boost::add_vertex(vp, *this);
    return std::pair<boost::shared_ptr<ClusterGraph>, LocalVertex> (m_clusters[v] = boost::shared_ptr<ClusterGraph> (new ClusterGraph(sp_base::shared_from_this())), v);
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
inline boost::shared_ptr< ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop> >
ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>:: parent() 	{
    return boost::shared_ptr<ClusterGraph> (m_parent);
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
inline const boost::shared_ptr< ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop> >
ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::parent() const 	{
    return boost::shared_ptr<ClusterGraph> (m_parent);
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
bool ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::isRoot() const {
    return m_parent.expired();
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
boost::shared_ptr< ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop> >
ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::root()		{
    return isRoot() ? sp_base::shared_from_this() : parent()->root();
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
const boost::shared_ptr< ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop> >
ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::root() const    {
    return isRoot() ? sp_base::shared_from_this() : parent()->root();
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
std::pair<typename ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::cluster_iterator, typename ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::cluster_iterator>
ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::clusters() {
    return std::make_pair(m_clusters.begin(), m_clusters.end());
}

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
std::pair<typename ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::const_cluster_iterator, typename ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::const_cluster_iterator>
ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::clusters() const {
    return std::make_pair(m_clusters.begin(), m_clusters.end());
}

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
std::size_t
ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::numClusters() const {
    return m_clusters.size();
}

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
bool ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::isCluster(const LocalVertex v) const {
    return (m_clusters.find(v) != m_clusters.end());
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
boost::shared_ptr< ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop> >
ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::getVertexCluster(LocalVertex v) {
    if(isCluster(v))
        return m_clusters[v];

    //TODO:throw if not a cluster
    return boost::shared_ptr< ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop> >();//sp_base::shared_from_this();
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
LocalVertex	ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::getClusterVertex(boost::shared_ptr<ClusterGraph> g) {
    std::pair<cluster_iterator, cluster_iterator> it = clusters();

    for(; it.first != it.second; it.first++) {
        if((*it.first).second == g)
            return (*it.first).first;
    }

    throw details::cluster_error() <<  boost::errinfo_errno(12) << error_message("Cluster is not part of this graph");
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
template<typename Functor>
void ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::removeCluster(boost::shared_ptr<ClusterGraph> g, Functor& f) {
    removeCluster(getClusterVertex(g), f);
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
void ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::removeCluster(boost::shared_ptr<ClusterGraph> g) {
    placehoder p;
    removeCluster(getClusterVertex(g), p);
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
void ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::clearClusters() {
    m_clusters.clear();
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
template<typename Functor>
void ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::removeCluster(LocalVertex v, Functor& f) {

    typename ClusterMap::iterator it = m_clusters.find(v);

    if(it == m_clusters.end())
        throw details::cluster_error() <<  boost::errinfo_errno(11) << error_message("Cluster is not part of this graph");

    std::pair<LocalVertex, boost::shared_ptr<ClusterGraph> > res = *it;

    //apply functor to all vertices and edges in the subclusters
    f(res.second);
    res.second->remove_vertices(f, true);

    //remove from map, delete subcluster and remove vertex
    m_clusters.erase(v);
    boost::clear_vertex(v, *this);    //should not be needed, just to ensure it
    boost::remove_vertex(v, *this);
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
void ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::removeCluster(LocalVertex v) {
    placehoder p;
    removeCluster(v, p);
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
template<typename Functor>
void ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::remove_vertices(Functor& f, bool recursive) {

    std::pair<local_vertex_iterator, local_vertex_iterator>  vit = boost::vertices(*this);

    //we iterate forward before deleting to not invalidate our iterator
    while(vit.first != vit.second) {
        LocalVertex v = * (vit.first);
        vit.first++;

        if(!isCluster(v)) {
            //let the functor know we remove this vertex
            f(getGlobalVertex(v));
            //need to do this to allow the removal of all relevant edges to this vertex, even upstream
            removeVertex(v, f);
        }
    };

    if(recursive) {
        cluster_iterator cit;

        for(cit = m_clusters.begin(); cit != m_clusters.end(); cit++) {
            f((*cit).second);
            (*cit).second->remove_vertices(f, recursive);
        }
    }
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
fusion::vector<LocalVertex, GlobalVertex>
ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::addVertex() {

    vertex_bundle vp;
    fusion::at_c<0> (vp) = m_id->generate();
    LocalVertex v = boost::add_vertex(vp, *this);

    return fusion::make_vector(v, m_id->count());
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
fusion::vector<LocalVertex, GlobalVertex>
ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::addVertex(GlobalVertex gv) {

    std::pair<LocalVertex, bool> res = getLocalVertex(gv);

    if(!res.second) {
        vertex_bundle vp;
        fusion::at_c<0> (vp) = gv;
        LocalVertex v = boost::add_vertex(vp, *this);

        //ensure that we never create this id, as it is used now
        if(gv > m_id->count())
            m_id->setCount(gv);

        return fusion::make_vector(v, gv);
    };

    return fusion::make_vector(res.first, gv);
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
std::pair<typename ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::global_vertex_iterator, typename ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::global_vertex_iterator>
ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::globalVertices() {
    std::pair<local_vertex_iterator, local_vertex_iterator> res = boost::vertices(*this);
    global_vertex_iterator begin = boost::make_transform_iterator(res.first, global_vertex_extractor(*this));
    global_vertex_iterator end   = boost::make_transform_iterator(res.second, global_vertex_extractor(*this));

    return std::pair<global_vertex_iterator, global_vertex_iterator> (begin, end);
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
std::pair<LocalEdge, bool>
ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::edge(LocalVertex source, LocalVertex target) {
    return boost::edge(source, target, *this);
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
fusion::vector<LocalEdge, GlobalEdge, bool>
ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::addEdge(LocalVertex source, LocalVertex target) {

    //manual edge creation with cluster is not allowed
    if((source == target) || isCluster(source) || isCluster(target))
        return fusion::make_vector(LocalEdge(), GlobalEdge(), false);

    LocalEdge e;
    bool done;
    boost::tie(e, done) = boost::edge(source, target, *this);

    //if done=true the edge alredy existed
    if(!done)
        boost::tie(e, done) = boost::add_edge(source, target, *this);

    if(!done)
        return fusion::make_vector(LocalEdge(), GlobalEdge(), false);

    //init the bundle corecctly for new edge
    GlobalEdge global = { fusion::at_c<0> ((*this) [source]), fusion::at_c<0> ((*this) [target]), m_id->generate() };
    edge_bundle_single s;
    fusion::at_c<0> (s) = global;
    fusion::at_c<1> ((*this) [e]).push_back(s);

    return fusion::make_vector(e, global, true);
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
fusion::vector<LocalEdge, GlobalEdge, bool, bool>
ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::addEdge(GlobalVertex source, GlobalVertex target) {

    LocalVertex v1, v2;
    LocalEdge e;
    bool d1, d2, d3;
    boost::tie(v1, d1) = getContainingVertex(source);
    boost::tie(v2, d2) = getContainingVertex(target);

    //if one vertex is not accessible from here this function fails
    if(!(d1 && d2))
        return fusion::make_vector(LocalEdge(), GlobalEdge(), false, false);

    //if both vertices are in a subcluster this one must do the job as we cant access the local edge from here
    if(v1 == v2 && isCluster(v1)) {
        fusion::vector<LocalEdge, GlobalEdge, bool, bool> res = getVertexCluster(v1)->addEdge(source, target);
        fusion::at_c<3> (res) = false;
        return res;
    }

    //check if we already have that Local edge
    boost::tie(e, d3) = boost::edge(v1, v2, *this);

    if(!d3)
        boost::tie(e, d3) = boost::add_edge(v1, v2, *this);

    if(!d3)
        return fusion::make_vector(LocalEdge(), GlobalEdge(), false, false);

    //init the bundle corectly for new edge
    GlobalEdge global = { source, target, m_id->generate() };
    edge_bundle_single s;
    fusion::at_c<0> (s) = global;
    fusion::at_c<1> ((*this) [e]).push_back(s);

    return fusion::make_vector(e, global, true, true);

};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
fusion::vector<LocalEdge, GlobalEdge, bool, bool>
ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::addEdgeGlobal(GlobalVertex source, GlobalVertex target) {
    return addEdge(source, target);
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
std::pair<typename ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::global_edge_iterator, typename ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::global_edge_iterator>
ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::getGlobalEdges(LocalEdge e) {

    std::vector<edge_bundle_single>& vec = fusion::at_c<1> ((*this) [e]);
    global_edge_iterator begin = boost::make_transform_iterator(vec.begin(), global_extractor());
    global_edge_iterator end   = boost::make_transform_iterator(vec.end(), global_extractor());

    return std::pair<global_edge_iterator, global_edge_iterator> (begin, end);
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
int ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::getGlobalEdgeCount(LocalEdge e) {

    return fusion::at_c<1> ((*this) [e]).size();
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
std::pair<LocalEdge, bool>
ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::getLocalEdge(GlobalEdge e) {
    return getContainingEdge(e);
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
fusion::vector<LocalEdge, ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>*, bool>
ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::getLocalEdgeGraph(GlobalEdge e) {
    return getContainingEdgeGraph(e);
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
GlobalVertex
ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::getGlobalVertex(LocalVertex v) const {
    return fusion::at_c<0> ((*this) [v]);
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
std::pair<LocalVertex, bool>
ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::getLocalVertex(GlobalVertex vertex) {
    return getContainingVertex(vertex);
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
fusion::vector<LocalVertex, boost::shared_ptr< ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop> >, bool>
ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::getLocalVertexGraph(GlobalVertex v) {
    return getContainingVertexGraph(v);
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
template<typename Functor>
void ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::downstreamRemoveVertex(GlobalVertex v, Functor& f) {

    std::pair<LocalVertex, bool> res = getContainingVertex(v);

    //we don't throw, as this function gets invoked recursivly and it may happen that the
    //vertex to remove is only in the top layers, not the button ones
    if(!res.second)
        return;


    //iterate over every edge that connects to the global vertex or the cluster in which it is in
    std::vector<LocalEdge> re; //remove edges
    std::pair<local_out_edge_iterator,  local_out_edge_iterator> it = boost::out_edges(res.first, *this);

    for(; it.first != it.second; it.first++) {
        std::vector<edge_bundle_single>& vec = fusion::at_c<1> ((*this) [* (it.first)]);
        vec.erase(std::remove_if(vec.begin(), vec.end(), apply_remove_prediacte<Functor, ClusterGraph> (f, v)), vec.end());

        if(vec.empty())
            re.push_back(* (it.first));
    };

    std::for_each(re.begin(), re.end(), boost::bind(&ClusterGraph::simpleRemoveEdge, this, _1));

    //if we have the real vertex here and not only a containing cluster we can delete it
    if(!isCluster(res.first)) {
        boost::clear_vertex(res.first, *this);    //just to make sure, should be done already
        boost::remove_vertex(res.first, *this);
    };

    //lets go downstream
    for(cluster_iterator it = m_clusters.begin(); it != m_clusters.end(); it++)
        ((*it).second)->downstreamRemoveVertex(v, f);
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
void ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::simpleRemoveEdge(LocalEdge e) {
    boost::remove_edge(e, *this);
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
template<typename Functor>
void ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::removeVertex(LocalVertex id, Functor& f) {
    //it is important to delete the global vertex, not the only local one as it's possible that
    //we are in a subcluster and there are connections to the global vertex in the parent. They
    //need to be deleted too.
    removeVertex(getGlobalVertex(id), f);
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
void ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::removeVertex(LocalVertex id) {
    placehoder p;
    removeVertex(getGlobalVertex(id), p);
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
template<typename Functor>
void ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::removeVertex(GlobalVertex id, Functor& f) {
    root()->downstreamRemoveVertex(id, f);
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
void ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::removeVertex(GlobalVertex id) {
    placehoder p;
    removeVertex(id, p);
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
void ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::removeEdge(GlobalEdge id) {
    fusion::vector<LocalEdge, ClusterGraph*, bool> res = getContainingEdgeGraph(id);

    if(!fusion::at_c<2> (res))
        return; //TODO:throw

    placehoder p;
    std::vector<edge_bundle_single>& vec = fusion::at_c<1> ((*fusion::at_c<1> (res)) [fusion::at_c<0> (res)]);
    vec.erase(std::remove_if(vec.begin(), vec.end(), apply_remove_prediacte<placehoder, ClusterGraph> (p, id)), vec.end());

    if(vec.empty())
        boost::remove_edge(fusion::at_c<0> (res), *fusion::at_c<1> (res));
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
template<typename Functor>
void ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::removeEdge(LocalEdge id, Functor& f) {

    std::vector<edge_bundle_single>& vec = fusion::at_c<1> ((*this) [id]);
    std::for_each(vec.begin(), vec.end(), boost::bind<void> (boost::ref(apply_remove_prediacte<placehoder, ClusterGraph> (f, -1)), _1));
    boost::remove_edge(id, *this);
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
template<typename property, typename key>
const typename property::type&
ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::getProperty(key k) {
    return apply_to_bundle(k, get_property_helper<property, ClusterGraph>());
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
template<typename property, typename key>
void ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::setProperty(key k, const typename property::type& val) {
    apply_to_bundle(k, set_property_helper<property, ClusterGraph>(val));
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
template<typename property, typename key>
bool ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::isPropertyChanged(key k) {
    return apply_to_bundle(k, change_property_helper<property, ClusterGraph>());
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
template<typename property, typename key>
void ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::acknowledgePropertyChange(key k) {
    apply_to_bundle(k, acknowledge_property_helper<property, ClusterGraph>());
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
template<typename key>
bool ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::hasPropertyChanges(key k) {
    return apply_to_bundle(k, changes_property_helper<ClusterGraph, key>());
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
template<typename key>
void ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::acknowledgePropertyChanges(key k) {
    apply_to_bundle(k, ack_changes_property_helper<ClusterGraph, key>());
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
void ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::initIndexMaps() {

    //just iterate over all edges and vertices and give them all a unique index
    std::pair<local_vertex_iterator, local_vertex_iterator>  vit = boost::vertices(*this);

    for(int c = 0; vit.first != vit.second; vit.first++, c++)
        setProperty<vertex_index>(*vit.first, c);

    std::pair<local_edge_iterator, local_edge_iterator>  eit = boost::edges(*this);

    for(int c = 0; eit.first != eit.second; eit.first++, c++)
        setProperty<edge_index>(*eit.first, c);
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
LocalVertex
ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::moveToSubcluster(LocalVertex v, boost::shared_ptr<ClusterGraph> cg) {

    LocalVertex cv = getClusterVertex(cg);
    return moveToSubcluster(v, cv, cg);
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
LocalVertex
ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::moveToSubcluster(LocalVertex v, LocalVertex Cluster) {

    boost::shared_ptr<ClusterGraph> cg = getVertexCluster(Cluster);
    return moveToSubcluster(v, Cluster, cg);
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
LocalVertex
ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::moveToSubcluster(LocalVertex v, LocalVertex Cluster, boost::shared_ptr<ClusterGraph> cg) {

    std::pair<local_out_edge_iterator, local_out_edge_iterator> it =  boost::out_edges(v, *this);

    /* add the later removed edges to the coressponding existing edges
     * (or create new edges between adjacent vertices of moved vertex and cluster).
     * also get the edge between cluster and vertex while iterating */
    for(; it.first != it.second; it.first++) {

        LocalVertex target = boost::target(*it.first, *this);

        if(target != Cluster) {

            //get or create the edge between the old edge target and the cluster
            LocalEdge e;
            bool done;
            boost::tie(e, done) = boost::edge(target, Cluster, *this);

            if(!done)
                boost::tie(e, done) = boost::add_edge(target, Cluster, *this);

            //if(!done) TODO: throw

            std::vector<edge_bundle_single>& ep = fusion::at_c<1> ((*this) [*it.first]);
            std::vector<edge_bundle_single>& nep = fusion::at_c<1> ((*this) [e]);
            nep.insert(nep.end(), ep.begin(), ep.end());
        }
    }

    /* Create new Vertex in Cluster and map the edge to vertices and clusters in the cluster
    * if a connection existed */
    LocalVertex nv = boost::add_vertex((*this) [v], *cg);

    //resort cluster parentship if needed
    if(isCluster(v)) {

        cg->m_clusters[nv] = m_clusters[v];
        cg->m_clusters[nv]->m_parent = cg;
        m_clusters.erase(v);
    }

    std::pair<LocalEdge, bool> moveedge = boost::edge(v, Cluster, *this);

    if(moveedge.second) {
        std::vector<edge_bundle_single>& vec = fusion::at_c<1> ((*this) [moveedge.first]);

        for(edge_single_iterator i = vec.begin(); i != vec.end(); i++) {

            //get the global vertex to which the global edge points and find the local vertex holding this
            //global one
            GlobalEdge global = global_extractor()(*i);
            GlobalVertex target;
            //bit cumbersome to support moving clusters
            target = (cg->getContainingVertex(global.source).first == nv) ? global.target : global.source;
            std::pair<LocalVertex, bool> res = cg->getContainingVertex(target);
            //if(!res.second) TODO: throw

            //get or create the edge between the new vertex and the target
            LocalEdge e;
            bool done;
            boost::tie(e, done) = boost::edge(nv, res.first, *cg);

            if(!done)
                boost::tie(e, done) = boost::add_edge(nv, res.first, *cg);

            //if(!done) TODO: throw

            //push the global edge to the local edge
            fusion::at_c<1> ((*cg) [e]).push_back(*i);
        };
    }

    //all global edges concerning the move vertex are processed and it is moved to the subcluster,
    //lets destroy it in the local cluster
    boost::clear_vertex(v, *this);
    boost::remove_vertex(v, *this);

    return nv;
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
LocalVertex
ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::moveToParent(LocalVertex v) {

    //if(isRoot()) TODO:throw

    //create new vertex
    vertex_bundle& vb = (*this) [v];
    LocalVertex nv = boost::add_vertex(vb, *parent());

    //regrouping if needed
    if(isCluster(v)) {
        parent()->m_clusters[nv] = m_clusters[v];
        parent()->m_clusters[nv]->m_parent = m_parent;
        m_clusters.erase(v);
    }

    GlobalVertex gv = fusion::at_c<0> (vb);

    //get all out_edges of this cluster in the parentcluster (because only they can hold relevant global_Edgs)
    std::vector<LocalEdge> edge_vec;
    LocalVertex this_v = parent()->getClusterVertex(sp_base::shared_from_this());
    std::pair<local_out_edge_iterator, local_out_edge_iterator> it = boost::out_edges(this_v, *parent());

    for(; it.first != it.second; it.first++) {
        //iterate all global edges and find relevant ones
        std::vector<edge_bundle_single>& vec = fusion::at_c<1> ((*parent()) [*it.first]);
        edge_single_iterator i = vec.begin();

        while(i != vec.end()) {

            GlobalEdge global = global_extractor()(*i);
            GlobalVertex target;

            //a bit cumbersome to allow cluster moving
            if(parent()->getContainingVertex(global.source).first == nv)
                target = global.target;
            else if(parent()->getContainingVertex(global.target).first == nv)
                target = global.source;
            else {
                i++;
                continue;
            }

            std::pair<LocalVertex, bool> res = parent()->getContainingVertex(target);

            //get or create the edge between the new vertex and the target
            LocalEdge e;
            bool done;
            boost::tie(e, done) = boost::edge(nv, res.first, *parent());

            if(!done)
                boost::tie(e, done) = boost::add_edge(nv, res.first, *parent());

            //if(!done) TODO: throw

            //push the global edge bundle to the new local edge and erase it in the old
            fusion::at_c<1> ((*parent()) [e]).push_back(*i);
            i = vec.erase(i);
        }

        //see if we should destroy this edge (no global edges remain in local one)
        if(vec.empty())
            edge_vec.push_back(*it.first);
    }

    //create a edge between new vertex and this cluster and add all global edges from within this cluster
    it = boost::out_edges(v, *this);
    LocalEdge e;

    if(it.first != it.second)
        e = boost::add_edge(nv, this_v, *parent()).first;

    for(; it.first != it.second; it.first++) {
        std::vector<edge_bundle_single>& ep = fusion::at_c<1> ((*this) [*it.first]);
        std::vector<edge_bundle_single>& nep = fusion::at_c<1> ((*parent()) [e]);
        nep.insert(nep.end(), ep.begin(), ep.end());
    }

    //all global edges concerning the move vertex are processed and it is moved to the parent,
    //lets destroy it in the local cluster
    boost::clear_vertex(v, *this);
    boost::remove_vertex(v, *this);

    //it's possible that some local edges in the parent are empty now, let's destroy them
    for(std::vector<LocalEdge>::iterator it = edge_vec.begin(); it != edge_vec.end(); it++)
        boost::remove_edge(*it, *parent());

    return nv;
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
std::pair<LocalVertex, bool>
ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::getContainingVertex(GlobalVertex id, bool recursive) {

    //check all vertices if they are the id
    std::pair<local_vertex_iterator, local_vertex_iterator>  it = boost::vertices(*this);

    for(; it.first != it.second; it.first++) {
        if(id == fusion::at_c<0> ((*this) [*it.first]))
            return std::make_pair(*it.first, true);
    }

    //check all clusters if they have the id
    if(recursive) {
        for(cluster_iterator it = m_clusters.begin(); it != m_clusters.end(); it++) {
            std::pair<LocalVertex, bool> res = ((*it).second)->getContainingVertex(id);

            if(res.second)
                return std::make_pair((*it).first, true);
        }
    }

    return std::make_pair((LocalVertex) NULL, false);
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
fusion::vector<LocalVertex, boost::shared_ptr< ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop> >, bool>
ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::getContainingVertexGraph(GlobalVertex id) {

    LocalVertex v;
    bool done;
    boost::tie(v, done) = getContainingVertex(id);

    if(!done)
        return fusion::make_vector(LocalVertex(), boost::shared_ptr<ClusterGraph>(), false);

    if(isCluster(v) && (getGlobalVertex(v) != id))
        return m_clusters[v]->getContainingVertexGraph(id);
    else
        return fusion::make_vector(v, sp_base::shared_from_this(), true);
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
std::pair<LocalEdge, bool>
ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::getContainingEdge(GlobalEdge id) {

    LocalVertex v1, v2;
    bool d1, d2;
    boost::tie(v1, d1) = getContainingVertex(id.source, true);
    boost::tie(v2, d2) = getContainingVertex(id.target, true);

    if(!((d1 && d2) && (v1 != v2)))
        return std::make_pair(LocalEdge(), false);

    return boost::edge(v1, v2, *this);
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
fusion::vector<LocalEdge, ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>*, bool>
ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::getContainingEdgeGraph(GlobalEdge id) {

    LocalVertex v1, v2;
    bool d1, d2;
    boost::tie(v1, d1) = getContainingVertex(id.source, true);
    boost::tie(v2, d2) = getContainingVertex(id.target, true);

    if(!(d1 && d2))
        return fusion::make_vector(LocalEdge(), (ClusterGraph*) NULL, false);

    if(v1 == v2)
        return m_clusters[v1]->getContainingEdgeGraph(id);

    return fusion::make_vector(boost::edge(v1, v2, *this).first, this, true);
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
template<typename functor>
typename functor::result_type
ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::apply_to_bundle(LocalVertex k, functor f) {
    return f((*this) [k]);
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
template<typename functor>
typename functor::result_type
ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::apply_to_bundle(LocalEdge k, functor f) {
    return f((*this) [k]);
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
template<typename functor>
typename functor::result_type
ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::apply_to_bundle(GlobalVertex k, functor f) {

    //check all vertices if they are the id
    std::pair<local_vertex_iterator, local_vertex_iterator>  it = boost::vertices(*this);

    for(; it.first != it.second; it.first++) {
        vertex_bundle& p = (*this) [*it.first];

        if(k == fusion::at_c<0> (p))
            return f(p);
    }

    //check all clusters if they have the object
    fusion::vector<LocalVertex, boost::shared_ptr<ClusterGraph>, bool> res = getContainingVertexGraph(k);

    if(!fusion::at_c<2> (res)) {
        //TODO: Throw (propeties return reference, but cant init a reference temporarily)
    }

    return fusion::at_c<1> (res)->template apply_to_bundle<functor> (k, f);
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
template<typename functor>
typename functor::result_type
ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::apply_to_bundle(GlobalEdge k, functor f) {

    LocalVertex v1, v2;
    bool d1, d2;
    boost::tie(v1, d1) = getContainingVertex(k.source);
    boost::tie(v2, d2) = getContainingVertex(k.target);

    if(!(d1 && d2)) {
        //TODO:Throw
    }

    if((v1 == v2) && isCluster(v1))
        return m_clusters[v1]->apply_to_bundle(k, f);
    else {
        LocalEdge e;
        bool done;
        boost::tie(e, done) = boost::edge(v1, v2, *this);
        //if(!done) TODO: throw, as there has to be a edge!

        //search the global one in the local edge
        std::vector< edge_bundle_single >& vec = fusion::at_c<1>((*this) [e]);
        edge_single_iterator it;

        for(it = vec.begin(); it != vec.end(); it++) {
            if(fusion::at_c<0>(*it) == k)
                return f(*it);
        };

        //TODO: throw, as there has to be a global edge and the following return statement is illegal
        return f(*it);
    };


};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
template<typename Predicate, typename Iterator>
std::pair<boost::filter_iterator<Predicate, Iterator>, boost::filter_iterator<Predicate, Iterator> >
ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::filterRange(std::pair<Iterator, Iterator> range) {

    Predicate p = Predicate(sp_base::shared_from_this());
    return std::make_pair(boost::make_filter_iterator(p, range.first, range.second),
                          boost::make_filter_iterator(p, range.second, range.second));
};

} //namespace dcm


#endif // CLUSTERGRAPH_HPP





