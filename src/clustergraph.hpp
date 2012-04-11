/*
    <one line to give the program's name and a brief idea of what it does.>
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

#ifndef CLUSTERGRAPH_HPP
#define CLUSTERGRAPH_HPP

#include <map>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>

#include <boost/mpl/transform.hpp>
#include <boost/mpl/find.hpp>
#include <boost/mpl/and.hpp>

#include <boost/utility/enable_if.hpp>

#include <boost/fusion/container/vector.hpp>
#include <boost/fusion/sequence.hpp>
#include <boost/fusion/algorithm.hpp>
#include <boost/fusion/include/make_vector.hpp>

namespace mpl = boost::mpl;
namespace fusion = boost::fusion;

namespace dcm {

typedef boost::adjacency_list_traits<boost::slistS, boost::slistS, boost::undirectedS> list_traits;
typedef typename list_traits::vertex_descriptor LocalVertex;
typedef typename list_traits::edge_descriptor LocalEdge;

typedef int universalID;

typedef universalID GlobalVertex;
typedef std::pair<universalID, universalID> GlobalEdge;

struct IDgen {
    universalID* counter;

    IDgen() {
        counter = new universalID(10);
    };
    universalID generate() {
       return ++(*counter);
    };
    universalID count() {
        return (*counter);
    };
};
typedef boost::shared_ptr<IDgen> IDpointer;

struct clear_ptr {
    template<typename T>
    void operator()(T t) const {
        t.reset();
    };
};

template<typename T>
struct spv {
    typedef typename mpl::transform<T, boost::shared_ptr<mpl::_1> >::type type;
};


template< typename edge_prop, typename vertex_prop, typename objects>
class ClusterGraph : boost::adjacency_list< boost::slistS, boost::slistS,
            boost::undirectedS,
            fusion::vector<vertex_prop, typename spv<objects>::type, GlobalVertex >,
            std::vector< fusion::vector<edge_prop, typename spv<objects>::type,
            GlobalEdge > > >	{

    typedef fusion::vector<edge_prop, typename spv<objects>::type, GlobalEdge > edge_bundle_single;
    typedef std::vector< edge_bundle_single > edge_bundle;
    typedef fusion::vector<vertex_prop, typename spv<objects>::type, GlobalVertex > vertex_bundle;
    typedef boost::adjacency_list< boost::slistS, boost::slistS,
    boost::undirectedS, vertex_bundle, edge_bundle > Graph;

    typedef typename boost::graph_traits<Graph>::vertex_iterator   vertex_iterator;
    typedef typename boost::graph_traits<Graph>::edge_iterator     edge_iterator;
    typedef typename boost::graph_traits<Graph>::out_edge_iterator out_edge_iterator;

    typedef std::map<LocalVertex,ClusterGraph*> ClusterMap;

public:
    typedef typename ClusterMap::iterator 	cluster_iterator;

public:

    ClusterGraph ( ClusterGraph* g = 0 ) : m_parent ( g ), m_id(new IDgen) {

        if (g) m_id = g->m_id;
    };

    ~ClusterGraph() {

        if ( !m_parent ) {
            for ( typename ClusterMap::iterator i = m_clusters.begin(); i != m_clusters.end(); ++i )  {
                delete ( *i ).second;
            }
        }
        //TODO: if parent exists all vertices have to be transfered to it;
    };

    template<typename T>
    bool operator==(const T &other) const {
        return this == &other;
    };

    template<typename T>
    bool operator!=(const T &other) const {
        return !(this == &other);
    };

    
    /* *******************************************************
     * Subclustering
     * *******************************************************/

    ClusterGraph& 	createCluster() {
        LocalVertex v = boost::add_vertex ( *this );
        return * ( m_clusters[v] = new ClusterGraph ( this ) );
    };

    ClusterGraph& 	parent() 	{
        return *m_parent;
    };
    const ClusterGraph& parent() const 	{
        return *m_parent;
    };
    bool 		isRoot() const {
        return m_parent == 0;
    };
    ClusterGraph& 	root()		{
        return isRoot() ? *this : m_parent->root();
    };
    const ClusterGraph& root() const    {
        return isRoot() ? *this : m_parent->root();
    };

    std::pair<cluster_iterator, cluster_iterator> clusters() {
        return std::make_pair ( m_clusters.begin(), m_clusters.end() );
    }
    std::size_t numClusters() const {
        return m_clusters.size();
    }

    bool isCluster ( LocalVertex v ) {
        return ( m_clusters.find ( v ) != m_clusters.end() );
    };

    ClusterGraph& 	getVertexCluster ( LocalVertex v ) {
        if ( isCluster ( v ) )
            return *m_clusters[v];
    };

    LocalVertex		getClusterVertex ( ClusterGraph& g ) {
        std::pair<cluster_iterator, cluster_iterator> it = g.clusters();
        for ( ;it.first!=it.second; ++it.first ) {
            if ( ( *it.first ).second == &g )
                return ( *it.first ).first;
        }
    };

    bool vertexToCluster ( LocalVertex v, ClusterGraph& cluster ) {

        LocalVertex cv = getClusterVertex ( cluster );
        return vertexToCluster ( v, cv, cluster );
    }

    bool vertexToCluster ( LocalVertex v, LocalVertex cluster ) {

        ClusterGraph& cl = getVertexCluster ( cluster );
        return vertexToCluster ( v, cluster, cl );
    }

    bool vertexToCluster ( LocalVertex v, LocalVertex Cluster, ClusterGraph& cg ) {

        if (isCluster(v)) return false;

        std::pair<out_edge_iterator, out_edge_iterator> it =  boost::out_edges(v, *this);

        /* add the later removed edges to the coressponding existing edges
         * (or create new edges between adjacent vertices of moved vertex and cluster).
         * also get the edge between cluster and vertex while iterating */
        std::pair<edge_bundle&,bool> clusterp;
        clusterp.second = false;
        for (;it.first!=it.second;it.first++) {

            LocalVertex target = boost::target(*it.first, *this);
            if (target != Cluster) {

                std::pair<LocalEdge, bool> res = boost::add_edge(target, Cluster, *this);
                edge_bundle& ep = (*this)[*it.first];
                edge_bundle& nep = (*this)[res.first];
                nep.insert(nep.end(), ep.begin(), ep.end());
            }
            else clusterp = std::make_pair( (*this)[*it.first], true );
        }

        /* Create new Vertex in Cluster and map the edge to vertices and clusters in the cluster
        * if a connection existed */
        LocalVertex nv = boost::add_vertex((*this)[v], cg);
        if (clusterp.second) {
            for (typename edge_bundle::iterator i = clusterp.begin(); i != clusterp.end(); i++) {

                std::pair<universalID, universalID> idpair = fusion::at_c<2>(*i);
                universalID target = (idpair.first == fusion::at_c<2>(cg[nv])) ? idpair.second : idpair.first;
                std::pair<LocalVertex, bool> res = getContainingVertex(target);
                LocalEdge ne = boost::add_edge(nv, res.first, cg);
                cg[ne].push_back(*i);
            };
        }
    }

     /* *******************************************************
     * Creation Handling
     * *******************************************************/

     /**
      * @brief Add a vertex to the local cluster
      *
      * @return fusion:vector< LocalVertex, GlobalVertex > with the local and global vertex descriptor
      **/
     fusion::vector<LocalVertex, GlobalVertex> addVertex() {

        vertex_bundle vp;
        fusion::at_c<2>(vp) = m_id->generate();
        LocalVertex v= boost::add_vertex(vp, *this);
        return fusion::make_vector(v, m_id->count());
    };

    
    /**
     * @brief Add a edge between two vertices, defined by local descriptors.
     * 
     * Add an edge that connects the two vertices and in the local clustergraph and assign the GlobalEdge to it. The
     * LocalVertex parameters should not represent a cluster which would result in the functions failure. If ther's
     * already a edge between the vertices the existing local and global descriptors are returnd. This counts as 
     * successful creation and will therefore not be indicated. Failure will be recocnisable by a false value in the 
     * returned type sequence.
     *
     * @param source The first vertex the edge should connect
     * @param target The second vertex the edge should connect
     * @return fusion::vector<LocalEdge, GlobalEdge, success> with the local and global descriptors of the edge and an bool
     * value indicationg the successful creation.
     **/
    fusion::vector<LocalEdge, GlobalEdge, bool> addEdge(LocalVertex source, LocalVertex target) {

	//manual edge creation with cluster is not allowed
	if( (source==target) || isCluster(source) || isCluster(target) ) 
	  return fusion::make_vector(LocalEdge(), GlobalEdge(), false);
	
        LocalEdge e;
        bool done;
        boost::tie(e,done) = boost::add_edge(source, target, *this);
	
	//if done=false the edge alredy existed
	if(!done) {
	  edge_bundle& vec = (*this)[e];
	  //if a non-cluster edge has more than one property attached something has gone terribly wrong
	  return fusion::make_vector(e, fusion::at_c<2>(vec.front()), vec.size()==1);
	}
	
	//init the bundle corecctly for new edge
	GlobalEdge global = std::make_pair(fusion::at_c<2>((*this)[source]),fusion::at_c<2>((*this)[target]));
	edge_bundle_single s;
	fusion::at_c<2>(s) = global;
	(*this)[e].push_back(s);
        return fusion::make_vector(e, global, true);
    };

    /**
     * @brief Add a edge between two vertices, defined by global descriptors.
     * 
     * Adds an edge between vertices which are not nesseccarily in this local cluster and have therefore to be 
     * identified with global descriptors. The only condition for source and tearget vertex is that both must be
     * in the local cluster or any of its child clusters. If thats not the case, the function will fail. On success
     * a new GlobalEdge will be created, but not neccessarily a local one. If the vertices are in different cluster
     * which are already connected the global edge will be added to this connecting local edge. Thats the one returned
     * in the seqence. Note that it's possible that the local edge belongs to another subcluster and therefore can't be
     * used in the local cluster. This case is indicated by the scope return value. 
     * In the case of an alredy existing global edge this one will be returned. Note that parent globaledges are not 
     * checked which may lead to existane of double edges, which would be a serious error.
     * 
     *
     * @param source The first vertex the edge should connect
     * @param target The second vertex the edge should connect
     * @return fusion:vector< LocalEdge, GlobalEdge, success, scope > with the new global edge descriptor and the local
     * one where it was added. Success indicates if the function was successful and scope shows the validy of the local
     * descriptor in this cluster (true means the edge is in this cluster).
     **/
    
    fusion::vector<LocalEdge, GlobalEdge, bool, bool> addEdge(GlobalVertex source, GlobalVertex target) {

        LocalVertex v1,v2;
        LocalEdge e;
        bool d1,d2,d3;
        boost::tie(v1,d1) = getContainingVertex(source);
        boost::tie(v2,d2) = getContainingVertex(target);
	
	//if one vertex is not accessible from here this function fails
	if (!(d1&&d2)) return fusion::make_vector(LocalEdge(), GlobalEdge(), false, false);
	
	//if both vertices are in a subcluster this one must do the job as we cant access the local edge from here
	if(v1==v2 && isCluster(v1)) {
	  fusion::vector<LocalEdge, GlobalEdge, bool, bool> res = getVertexCluster(v1)->addEdge(source, target);
	  fusion::at_c<3>(res)=false;
	  return res;
	}
	
        boost::tie(e,d3) = boost::add_edge(v1, v2, *this);
        //init the bundle corectly for new edge
	GlobalEdge global = std::make_pair(source, target);
	edge_bundle_single s;
	fusion::at_c<2>(s) = global;
	(*this)[e].push_back(s);
        return fusion::make_vector(e, global, true, true);
	
	//TODO: check if global edge already exists
    };

    /**
     * @brief Get the local edges global descriptor
     * 
     * Gives the global descriptor of the local edge. Note that local edges can connect to atleast one cluster,
     * in which case they can have multiple global descriptors. In this case the first one will be returned. 
     *
     * @param e the local edge for which the global descriptor is wanted
     * @return GlobalEdge
     **/
    GlobalEdge getGlobalEdge(LocalEdge e) {
	edge_bundle vec = (*this)[e];
	if(vec.empty()) return GlobalEdge();
        return fusion::at_c<2>(*((*this)[e]).begin());
    };
    /**
     * @brief Get the local edge which holds the specified global edge.
     *
     * Note that GlobalEdge must be in a local edge of this cluster, means the connected vertices must be in this 
     * ore one of it's subclusters. Also if the containing LocalEdge is no in this cluster, but in one of it's 
     * subclusters, the function fails and the returned edge is invalid.
     * 
     * @param e GlobalEdge for which the containing local one is wanted
     * @return std:pair< LocalEdge, bool > whith the containing LocalEdge and a bool indicator if function was successful.
     **/    
    std::pair<LocalEdge, bool> getLocalEdge(GlobalEdge e) {
        return getContainingEdge(e, true);
    };
    
    /**
     * @brief Get the GlobalVertex assiociated with this local one. 
     *
     * @param e LocalVertex
     * @return GlobalVertex
     **/
    
    GlobalVertex getGlobalVertex(LocalVertex e) {
        return fusion::at_c<2>((*this)[e]);
    };
    /**
     * @brief Get the LocalVertex which corresponds to the golab one
     * 
     * The GlobalVertex has to be in this cluster or any of it's subclusters. If its in a subcluster, the returned 
     * LocalVertex will represent this cluster. If the GlobalVertex is not in this clusters scope the function fails.
     *
     * @param GlobalVertex
     * @return std::pair< LocalVertex, bool > The LocalVertex containing the global one and an success indicator
     **/    
    std::pair<LocalVertex,bool> getLocalVertex(GlobalVertex e) {
        return getContainingVertex(e);
    };


    /* *******************************************************
     * Object Handling
     * *******************************************************/
    template<typename Obj, typename key>
    struct get_obj_helper {

        get_obj_helper(key k) : m_key(k) {};

        typedef boost::shared_ptr<Obj> result_type;
        typedef typename mpl::find<objects, Obj>::type iterator;
        BOOST_MPL_ASSERT(( mpl::not_<boost::is_same<iterator, typename mpl::end<objects>::type > > ));

        template<typename prop>
        typename boost::enable_if<boost::is_same<prop, typename boost::vertex_bundle_type<Graph>::type>,
        result_type>::type operator()(prop& p) {
            return fusion::at_c<iterator::pos::value>(fusion::at_c<1>(p));
        }

        template<typename prop>
        typename boost::enable_if<mpl::and_<boost::is_same<prop, typename boost::edge_bundle_type<Graph>::type>,
        boost::is_same<key, GlobalEdge> >, result_type>::type operator()(prop& p) {

            for (typename prop::iterator it=p.begin(); it != p.end(); it++) {
                if ( fusion::at_c<2>(*it) == m_key )
                    return fusion::at_c<iterator::pos::value>(fusion::at_c<1>(*it));
            }
        }

        template<typename prop>
        typename boost::enable_if<mpl::and_<boost::is_same<prop, typename boost::edge_bundle_type<Graph>::type>,
        boost::is_same<key, LocalEdge> >, result_type>::type operator()(prop& p) {
            return fusion::at_c<iterator::pos::value>(fusion::at_c<1>(*p.begin()));
        }

        key m_key;
    };

    template<typename Obj, typename key>
    boost::shared_ptr<Obj> getObject(key k) {
        return apply_to_bundle(k, get_obj_helper<Obj, key>(k));
    };

    template<typename Obj, typename key>
    void setObject(key k, boost::shared_ptr<Obj> val) {
        //apply_to_bundle(k, get_obj_helper<Obj>()) = val;
    };



    /* *******************************************************
     * Property Handling
     * *******************************************************/
    template<typename prop, typename key>
    struct get_prop_helper {

        get_prop_helper(key k) : m_key(k) {};

	typedef typename prop::type& result_type;
        typedef typename mpl::find<vertex_prop, prop>::type vertex_iterator;
        typedef typename mpl::find<edge_prop, prop>::type edge_iterator;
        typedef typename mpl::if_<boost::is_same<vertex_iterator,typename mpl::end<vertex_prop>::type >,
        edge_iterator, vertex_iterator>::type iterator;
	BOOST_MPL_ASSERT(( mpl::not_<boost::is_same<iterator, typename mpl::end<edge_prop>::type > > ));

        template<typename bundle>
        typename boost::enable_if<boost::is_same<prop, typename boost::vertex_bundle_type<Graph>::type>,
        result_type>::type operator()(prop& p) {
            return fusion::at_c<iterator::pos::value>(fusion::at_c<0>(p));
        }

        template<typename bundle>
        typename boost::enable_if<mpl::and_<boost::is_same<prop, typename boost::edge_bundle_type<Graph>::type>,
        boost::is_same<key, GlobalEdge> >, result_type>::type operator()(prop& p) {

            for (typename bundle::iterator it=p.begin(); it != p.end(); it++) {
                if ( fusion::at_c<2>(*it) == m_key )
                    return fusion::at_c<iterator::pos::value>(fusion::at_c<0>(*it));
            }
        }

        template<typename bundle>
        typename boost::enable_if<mpl::and_<boost::is_same<prop, typename boost::edge_bundle_type<Graph>::type>,
        boost::is_same<key, LocalEdge> >, result_type>::type operator()(prop& p) {
            return fusion::at_c<iterator::pos::value>(fusion::at_c<0>(*p.begin()));
        }

        key m_key;
    };

    template<typename prop, typename key>
    typename prop::type& getProperty(key k) {
        return apply_to_bundle(k, get_prop_helper<prop, key>(k));
    };

    template<typename prop, typename key>
    void setProperty(key k, prop val) {
        apply_to_bundle(k, get_prop_helper<prop, key>(k)) = val;
    };

protected:
    ClusterGraph* m_parent;
    ClusterMap	  m_clusters;
    IDpointer 	  m_id;

    std::pair<LocalVertex, bool> getContainingVertex(GlobalVertex id, bool recursive = true) {

        //check all vertices if they are the id
        std::pair<vertex_iterator, vertex_iterator>  it = boost::vertices(*this);
        for (;it.first != it.second; it.first++) {
            if ( id == fusion::at_c<2>((*this)[*it.first]) )
               return std::make_pair(*it.first, true);
        }

        //check all clusters if they have the id
        if (recursive) {
            for (cluster_iterator it = m_clusters.begin(); it != m_clusters.end(); it++) {
                std::pair<LocalVertex, bool> res =  ((*it).second)->getContainingVertex(id);
                if (res.second) return std::make_pair((*it).first, true);
            }
        }

        return std::make_pair((LocalVertex)NULL, false);
    };

    std::pair<LocalEdge, bool> getContainingEdge(GlobalEdge id, bool recusive = true) {

        LocalVertex v1,v2;
        bool d1,d2;
        boost::tie(v1,d1) = getContainingVertex(id.first, recusive);
        boost::tie(v2,d2) = getContainingVertex(id.second, recusive);

        if ( !((d1&&d2) && (v1!=v2)) )   return std::make_pair(LocalEdge(), false);

        return boost::edge(v1,v2,*this);
    };

    std::pair<LocalEdge, ClusterGraph*> getContainingEdgeGraph(GlobalEdge id) {

        LocalVertex v1,v2;
        bool d1,d2;
        boost::tie(v1,d1) = getContainingVertex(id.first);
        boost::tie(v2,d2) = getContainingVertex(id.second);

        if ( d1&&d2 ) {
            if ( (v1==v2) && isCluster(v1) ) return m_clusters[v1]->getContainingEdgeGraph(id);
            else return std::make_pair(boost::edge(v1,v2, *this).first, this);
        }

        return std::make_pair(NULL, NULL);
    };

    template<typename functor>
    typename functor::result_type apply_to_bundle(LocalVertex k, functor f) {
        f((*this)[k]);
    };
    template<typename functor>
    typename functor::result_type apply_to_bundle(LocalEdge k, functor f) {
        f((*this)[k]);
    };

    template<typename functor>
    typename functor::result_type apply_to_bundle(GlobalVertex k, functor f) {

        //check all vertices if they are the id
        std::pair<vertex_iterator, vertex_iterator>  it = boost::vertices(*this);
        for (;it.first != it.second; it.first++) {
            vertex_prop& p = (*this)[*it.first];
            if ( k == fusion::at_c<2>(p) )
                return f(p);
        }

        //check all clusters if they have the object
        for (cluster_iterator it = m_clusters.begin(); it != m_clusters.end(); it++) {
            return (*it.second)->apply_to_bundle<functor>(k, f);
        }

        return functor::result_type();
    };

    template<typename functor>
    typename functor::result_type apply_to_bundle(GlobalEdge k, functor f) {

        std::pair<LocalEdge, ClusterGraph*> res = getContainingEdgeGraph(k);
        if (!res.second) return functor::result_type();
        return f((*this)[k]);
    };
};

} //namespace solver

#endif // CLUSTERGRAPH_HPP
