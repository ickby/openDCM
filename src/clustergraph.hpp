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

#ifndef CLUSTERGRAPH_HPP
#define CLUSTERGRAPH_HPP

#include <map>
#include <functional>
#include <iostream>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>

#include <boost/mpl/transform.hpp>
#include <boost/mpl/find.hpp>
#include <boost/mpl/and.hpp>

#include <boost/utility/enable_if.hpp>
#include <boost/iterator/transform_iterator.hpp>

#include <boost/fusion/include/make_vector.hpp>
#include <boost/fusion/include/at_c.hpp>
#include <boost/fusion/include/at.hpp>
#include <boost/fusion/include/as_vector.hpp>
#include <boost/fusion/include/mpl.hpp>
#include <boost/fusion/container/vector.hpp>
#include <boost/fusion/sequence.hpp>
#include <boost/fusion/algorithm.hpp>

namespace mpl = boost::mpl;
namespace fusion = boost::fusion;

namespace dcm {

namespace details {

typedef boost::adjacency_list_traits<boost::slistS, boost::slistS, boost::undirectedS> list_traits;
typedef int universalID;

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
    void operator()(T& t) const {
        t.reset();
    };
};

template<typename T>
struct sps { //shared_ptr sequence
    typedef typename mpl::transform<T, boost::shared_ptr<mpl::_1> >::type spv;
    typedef typename fusion::result_of::as_vector<spv>::type type;
};

template<typename T>
struct prop_type {
    typedef typename T::type type;
};

template<typename T>
struct pts {
    typedef typename mpl::transform<T, prop_type<mpl::_1> >::type ptv;
    typedef typename fusion::result_of::as_vector< ptv >::type type;
};

}

typedef typename details::list_traits::vertex_descriptor 	LocalVertex;
typedef typename details::list_traits::edge_descriptor 		LocalEdge;
typedef details::universalID 					GlobalVertex;
typedef std::pair<details::universalID, details::universalID> 	GlobalEdge;


template< typename edge_prop, typename vertex_prop, typename objects>
class ClusterGraph : public boost::adjacency_list< boost::slistS, boost::slistS,
            boost::undirectedS,
            fusion::vector< typename details::pts<vertex_prop>::type,
            typename details::sps<objects>::type, GlobalVertex >,
            std::vector< fusion::vector<typename details::pts<edge_prop>::type,
            typename details::sps<objects>::type, GlobalEdge > > >	{

    typedef fusion::vector< typename details::pts<edge_prop>::type,
    typename details::sps<objects>::type, GlobalEdge > edge_bundle_single;
    typedef std::vector< edge_bundle_single > edge_bundle;
    typedef fusion::vector< typename details::pts<vertex_prop>::type,
    typename details::sps<objects>::type, GlobalVertex > vertex_bundle;
    typedef boost::adjacency_list< boost::slistS, boost::slistS,
    boost::undirectedS, vertex_bundle, edge_bundle > Graph;

    typedef typename boost::graph_traits<Graph>::vertex_iterator   local_vertex_iterator;
    typedef typename boost::graph_traits<Graph>::edge_iterator     local_edge_iterator;
    typedef typename boost::graph_traits<Graph>::out_edge_iterator local_out_edge_iterator;

    typedef std::map<LocalVertex,ClusterGraph*> ClusterMap;

private:
    struct global_extractor  {
        typedef GlobalEdge& result_type;
        template<typename T>
        result_type operator() (T& bundle) const {
            return fusion::at_c<2>(bundle);
        };
    };
    template<typename Obj>
    struct object_extractor  {

        typedef boost::shared_ptr<Obj> base_type;
        typedef base_type& result_type;
        typedef typename mpl::find<objects, Obj>::type iterator;
        typedef typename mpl::distance<typename mpl::begin<objects>::type, iterator>::type distance;
        BOOST_MPL_ASSERT(( mpl::not_<boost::is_same<iterator, typename mpl::end<objects>::type > > ));

        template<typename seq>
        result_type operator() (seq& bundle) const {
            return fusion::at<distance>(fusion::at< mpl::int_<1> >(bundle));
        };
    };

    template<typename prop>
    struct property_extractor  {

        typedef typename prop::type base_type;
        typedef base_type& result_type;
        typedef typename mpl::find<vertex_prop, prop>::type vertex_iterator;
        typedef typename mpl::find<edge_prop, prop>::type edge_iterator;
        typedef typename mpl::if_<boost::is_same<vertex_iterator,typename mpl::end<vertex_prop>::type >,
        edge_iterator, vertex_iterator>::type iterator;
        typedef typename mpl::if_<boost::is_same<vertex_iterator,typename mpl::end<vertex_prop>::type >,
        typename mpl::begin<edge_prop>::type, typename mpl::begin<vertex_prop>::type>::type begin;
        typedef typename mpl::distance<begin, iterator>::type distance;

        template<typename seq>
        result_type operator() (seq& bundle) const {
            return fusion::at<distance>(fusion::at< mpl::int_<0> >(bundle));
        };
    };

public:
    //iterators
    typedef boost::transform_iterator<global_extractor, typename edge_bundle::iterator> global_edge_iterator;

    template<typename Obj>
    struct object_iterator : public boost::transform_iterator<object_extractor<Obj>,typename edge_bundle::iterator> {
        object_iterator(typename edge_bundle::iterator it, object_extractor<Obj> f)
                : boost::transform_iterator<object_extractor<Obj>,typename edge_bundle::iterator>(it, f) {};
    };

    template<typename Prop>
    struct property_iterator : public boost::transform_iterator<property_extractor<Prop>,typename edge_bundle::iterator> {
        property_iterator(typename edge_bundle::iterator it, property_extractor<Prop> f)
                : boost::transform_iterator<property_extractor<Prop>,typename edge_bundle::iterator>(it, f) {};
    };


public:
    typedef typename ClusterMap::iterator 	cluster_iterator;

public:

    ClusterGraph ( ClusterGraph* g = 0 ) : m_parent ( g ), m_id(new details::IDgen) {

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

    //debug function to print the full typename
    template<typename T>
    void pretty (T t) {
        std::cout << __PRETTY_FUNCTION__ << std::endl;
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

    std::pair<ClusterGraph&, LocalVertex> createCluster() {
        LocalVertex v = boost::add_vertex ( *this );
        return std::pair<ClusterGraph&, LocalVertex>(* ( m_clusters[v] = new ClusterGraph ( this ) ), v);
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
        std::pair<cluster_iterator, cluster_iterator> it = clusters();
        for ( ;it.first!=it.second; it.first++ ) {
            if ( ( *it.first ).second == &g )
                return ( *it.first ).first;
        }
        return LocalVertex();
    };


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
     * @brief Returns the edge between the local vertices
     *
     * This function is the same as boost::edge(source, target, Graph) and only added for convienience.
     *
     * @param source LocalEdge as edge source
     * @param target LocalEdge as edge target
     * @return std::pair<LocalEdge, bool> with the local edge descriptor if existing. The bool value shows if the
     * edge exists or not
     **/
    std::pair<LocalEdge, bool> edge(LocalVertex source, LocalVertex target) {
        return boost::edge(source, target, *this);
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
        if ( (source==target) || isCluster(source) || isCluster(target) )
            return fusion::make_vector(LocalEdge(), GlobalEdge(), false);

        LocalEdge e;
        bool done;
        boost::tie(e,done) = boost::edge(source, target, *this);

        //if done=true the edge alredy existed
        if (done) {
            edge_bundle& vec = (*this)[e];
            //if a non-cluster edge has more than one property attached something has gone terribly wrong
            return fusion::make_vector(e, fusion::at_c<2>(vec.front()), vec.size()==1);
        }
        else boost::tie(e,done) = boost::add_edge(source, target, *this);

        if (!done) return fusion::make_vector(LocalEdge(), GlobalEdge(), false);

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
        if (v1==v2 && isCluster(v1)) {
            fusion::vector<LocalEdge, GlobalEdge, bool, bool> res = getVertexCluster(v1).addEdge(source, target);
            fusion::at_c<3>(res)=false;
            return res;
        }

        //check if we already have that Local and GlobalEdge
        boost::tie(e,d3) = boost::edge(v1,v2, *this);
        GlobalEdge global = std::make_pair(source, target);
        if (d3) {
            edge_bundle vec = (*this)[e];
            for (typename edge_bundle::iterator it=vec.begin(); it!=vec.end(); it++) {
                if (fusion::at_c<2>(*it) ==  global)
                    return fusion::make_vector(e, global, true, true);
            }
        }
        //local edge is not exsitant
        else boost::tie(e,d3) = boost::add_edge(v1, v2, *this);

        if (!d3) return fusion::make_vector(LocalEdge(), GlobalEdge(), false, false);

        //init the bundle corectly for new edge as the global edge is not existant
        edge_bundle_single s;
        fusion::at_c<2>(s) = global;
        (*this)[e].push_back(s);
        return fusion::make_vector(e, global, true, true);

    };

    /**
     * @brief Get an iterator to all the global edges hold by this local edge
     *
     * Local edges can hold multiple global ones, for example when they connect at least one cluster. Therefore a direct
     * LocalEdge - GlobalEdge mapping is not possible. Instead you can access all GlobalEdge's hold by this local one in
     * a normal iterating manner.
     *
     * @param e the local edge for which the global descriptors are wanted
     * @return std::pair<begin, end> with the global_edge_iterator's pointing to the vector<GlobalEdge>'s start
     * and end
     **/
    std::pair<global_edge_iterator, global_edge_iterator> getGlobalEdges(LocalEdge e) {

        edge_bundle& vec = (*this)[e];
        global_edge_iterator begin = boost::make_transform_iterator(vec.begin(), global_extractor());
        global_edge_iterator end   = boost::make_transform_iterator(vec.end(), global_extractor() );
        return std::pair<global_edge_iterator, global_edge_iterator>(begin, end);
    };


    /**
     * @brief Get the local edge which holds the specified global edge.
     *
     * Note that GlobalEdge must be in a local edge of this cluster, means the connected vertices must be in this
     * ore one of it's subclusters. Also if the containing LocalEdge is not in this cluster, but in one of it's
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

protected:
    //types needed to distinguish when objects need to be reset
    struct get : public boost::false_type {};
    struct set : public boost::true_type {};

    template<typename Type, typename Obj, typename key>
    struct obj_helper {

        typedef typename object_extractor<Obj>::result_type result_type;

        obj_helper(key k) : m_key(k) {};

        //used with vertex bundle type
        template<typename prop>
        typename boost::enable_if<boost::is_same<prop, typename boost::vertex_bundle_type<Graph>::type>,
        result_type >::type operator()(prop& p) {

            if (Type::value) fusion::for_each(fusion::at< mpl::int_<1> >(p), details::clear_ptr());
            return object_extractor<Obj>()(p);
        }

        //used with edge bundle type and global edge descriptor
        template<typename prop>
        typename boost::enable_if<mpl::and_<boost::is_same<prop, typename boost::edge_bundle_type<Graph>::type>,
        boost::is_same<key, GlobalEdge> >, result_type>::type operator()(prop& p) {

            //need to search the edge_bundle for the global descriptor
            for (typename prop::iterator it=p.begin(); it != p.end(); it++) {
                if ( fusion::at_c<2>(*it) == m_key ) {
                    if (Type::value) fusion::for_each(fusion::at< mpl::int_<1> >(*it), details::clear_ptr());
                    return object_extractor<Obj>()(*it);
                }
            }
        }

        //used with edge bundle type and local edge descriptor
        template<typename prop>
        typename boost::enable_if<mpl::and_<boost::is_same<prop, typename boost::edge_bundle_type<Graph>::type>,
        boost::is_same<key, LocalEdge> >, result_type>::type operator()(prop& p) {
            if (Type::value) fusion::for_each(fusion::at< mpl::int_<1> >(p.front()), details::clear_ptr());
            return object_extractor<Obj>()(p.front());
        }

        key m_key;
    };

public:

    /**
    * @brief Get the desired object at the specified vertex or edge
    *
    * This function allows to access the objects stored in the graph. If no object of the desired type
    * was set before, a empty shared_ptr will be returned. Accessing the object at a local edge is a special
    * case, as it can hold many global edges, each with it's own objetcs. Using a LocalEdge as key will
    * always return the object for the first GlobalEdge.
    *
    * @param local or global Vertex/Edge descriptor for which the object is desired
    * @return shared_ptr< Obj > the pointer to the desired object
    **/
    template<typename Obj, typename key>
    boost::shared_ptr<Obj> getObject(key k) {
        return apply_to_bundle(k, obj_helper<get, Obj, key>(k));
    };

    /**
     * @brief Set a object at the specified vertex or edge
     *
     * Sets the given value at the given key. Note that every entity can hold only one object, so setting
     * a new value resets all other objects which were set before. Setting the object at a local edge is a special
     * case, as it can hold many global edges, each with it's own objects. Using a LocalEdge as key will
     * always set the object for the first GlobalEdge.
     *
     * @param k local or global Vertex/Edge descriptor for which the object should be set
     * @param val the object which should be stored
     * @return void
     **/
    template<typename Obj, typename key>
    void setObject(key k, boost::shared_ptr<Obj> val) {
        apply_to_bundle(k, obj_helper<set, Obj, key>(k)) = val;
    };

    /**
     * @brief Get iterator range for all GlobalEdge objects hold by this local edge
     *
     * LocalEdge's can hold multiple global ones and the iterators can be used to access a specific object in
     * all global edges hold by this local edge.
     *
     * @param k the LocalEdge over which all Objects should be iterated.
     * @return pair< begin, end > the iterator rang from begin (first element) to end (first undefined element)
     **/
    template<typename Obj>
    std::pair< object_iterator<Obj>, object_iterator<Obj> > getObjects(LocalEdge k) {

        edge_bundle& vec = (*this)[k];
        object_iterator<Obj> begin(vec.begin(), object_extractor<Obj>());
        object_iterator<Obj> end(vec.end(), object_extractor<Obj>() );
        return std::pair< object_iterator<Obj>, object_iterator<Obj> >(begin, end);
    };




    /* *******************************************************
     * Property Handling
     * *******************************************************/

protected:

    template<typename prop, typename key>
    struct get_prop_helper {

        get_prop_helper(key k) : m_key(k) {};

        typedef typename prop::type base_type;
        typedef base_type& result_type;
        typedef typename mpl::find<vertex_prop, prop>::type vertex_iterator;
        typedef typename mpl::find<edge_prop, prop>::type edge_iterator;
        typedef typename mpl::if_<boost::is_same<vertex_iterator,typename mpl::end<vertex_prop>::type >,
        edge_iterator, vertex_iterator>::type iterator;
        BOOST_MPL_ASSERT(( mpl::not_<boost::is_same<iterator, typename mpl::end<edge_prop>::type > > ));

        template<typename bundle>
        typename boost::enable_if<boost::is_same<bundle, typename boost::vertex_bundle_type<Graph>::type>,
        result_type>::type operator()(bundle& p) {
            return property_extractor<prop>()(p);
        }

        template<typename bundle>
        typename boost::enable_if<mpl::and_<boost::is_same<bundle, typename boost::edge_bundle_type<Graph>::type>,
        boost::is_same<key, GlobalEdge> >, result_type>::type operator()(bundle& p) {

            for (typename bundle::iterator it=p.begin(); it != p.end(); it++) {
                if ( fusion::at_c<2>(*it) == m_key )
                    return property_extractor<prop>()(*it);
            }
        }

        template<typename bundle>
        typename boost::enable_if<mpl::and_<boost::is_same<bundle, typename boost::edge_bundle_type<Graph>::type>,
        boost::is_same<key, LocalEdge> >, result_type>::type operator()(bundle& p) {
            return property_extractor<prop>()(p.front());
        }

        key m_key;
    };

public:
    /**
    * @brief Get the desired property at the specified vertex or edge
    *
    * This function allows to access the properties stored in the graph. If no property of the desired type
    * was set before, a default construced will be returned. Accessing the property at a local edge is a special
    * case, as it can hold many global edges, each with it's own propertys. Using a LocalEdge as key will
    * always return the property for the first GlobalEdge.
    *
    * @param local or global Vertex/Edge descriptor for which the property is desired
    * @return property::type& the reference to the desired property
    **/
    template<typename property, typename key>
    typename property::type& getProperty(key k) {
        return apply_to_bundle(k, get_prop_helper<property, key>(k));
    };

    /**
     * @brief Set a property at the specified vertex or edge
     *
     * Sets the given value at the given key. Note that every entity can hold one of each property, as opposed
     * to objects. Setting the property at a local edge is a special case, as it can hold many global edges,
     * each with it's own propertys. Using a LocalEdge as key will always set the property for the first GlobalEdge.
     *
     * @param k local or global Vertex/Edge descriptor for which the property should be set
     * @param val the property value which should be stored
     * @return void
     **/
    template<typename property, typename key>
    void setProperty(key k, typename property::type val) {
        apply_to_bundle(k, get_prop_helper<property, key>(k)) = val;
    };

    /**
     * @brief Get iterator range for all GlobalEdge proeprties hold by this local edge
     *
     * LocalEdge's can hold multiple global ones and the iterators can be used to access a specific property in
     * all global edges hold by this local edge.
     *
     * @param k the LocalEdge over which all properties should be iterated.
     * @return pair< begin, end > the iterator rang from begin (first element) to end (first undefined element)
     **/
    template<typename property>
    std::pair< property_iterator<property>, property_iterator<property> > getProperties(LocalEdge k) {

        edge_bundle& vec = (*this)[k];
        property_iterator<property> begin(vec.begin(), property_extractor<property>());
        property_iterator<property> end(vec.end(), property_extractor<property>() );
        return std::pair< property_iterator<property>, property_iterator<property> >(begin, end);
    };



    /********************************************************
    * Vertex and Cluster moving
    * *******************************************************/

    /**
     * @brief Move a vertex to a subcluster
     *
     * Overloaded convinience function which fetches the local descriptor for the cluster reference and calls
     * the full parameter equivalent.
     *
     * @param v the LocalVertex to be moved
     * @param cg reference to the subcluster to which v should be moved
     * @return LocalVertex the local descriptor of the moved edge in the subcluster
     **/
    LocalVertex vertexToSubcluster ( LocalVertex v, ClusterGraph& cg ) {

        LocalVertex cv = getClusterVertex ( cg );
        return vertexToSubcluster ( v, cv, cg );
    };

    /**
     * @brief Move a vertex to a subcluster
     *
     * Overloaded convinience function which fetches the the cluster reference for the local descriptor and calls
     * the full parameter equivalent.
     *
     * @param v the LocalVertex to be moved
     * @param Cluster the local vertex descriptor representing the subcluster to which v should be moved
     * @return LocalVertex the local descriptor of the moved edge in the subcluster
     **/
    LocalVertex vertexToSubcluster ( LocalVertex v, LocalVertex Cluster ) {

        ClusterGraph& cg = getVertexCluster ( Cluster );
        return vertexToSubcluster ( v, Cluster, cg );
    };

    /**
     * @brief Move a vertex to a subcluster
     *
     * This function moves the LocalVertex to the subcluster and reconnects all other vertices and clusters. The
     * moved vertex will hold it's global descriptor but get a new local one assigned (the one returned). The same
     * stands for all edges which use the moved vertex: global descriptors stay the same, but they are moved to new
     * local edges.
     * The specified cluster has of course to be a valid subcluster.
     *
     * @param v the LocalVertex to be moved
     * @param Cluster the local vertex descriptor representing the subcluster to which v should be moved
     * @param cg reference to the subcluster to which v should be moved
     * @return LocalVertex the local descriptor of the moved edge in the subcluster
     **/
    LocalVertex vertexToSubcluster ( LocalVertex v, LocalVertex Cluster, ClusterGraph& cg ) {

        //if (isCluster(v)) return LocalVertex(); TODO: throw

        std::pair<local_out_edge_iterator, local_out_edge_iterator> it =  boost::out_edges(v, *this);

        /* add the later removed edges to the coressponding existing edges
         * (or create new edges between adjacent vertices of moved vertex and cluster).
         * also get the edge between cluster and vertex while iterating */
        for (;it.first!=it.second;it.first++) {

            LocalVertex target = boost::target(*it.first, *this);
            if (target != Cluster) {

                //get or create the edge between the old edge target and the cluster
                LocalEdge e;
                bool done;
                boost::tie(e,done) = boost::edge(target, Cluster, *this);
                if (!done ) boost::tie(e,done) = boost::add_edge(target, Cluster, *this);
                //if(!done) TODO: throw

                edge_bundle& ep = (*this)[*it.first];
                edge_bundle& nep = (*this)[e];
                nep.insert(nep.end(), ep.begin(), ep.end());
            }
        }

        /* Create new Vertex in Cluster and map the edge to vertices and clusters in the cluster
        * if a connection existed */
        LocalVertex nv = boost::add_vertex((*this)[v], cg);
        std::pair<LocalEdge, bool> moveedge = boost::edge(v, Cluster, *this);
        if (moveedge.second) {
            edge_bundle& vec = (*this)[moveedge.first];
            for (typename edge_bundle::iterator i = vec.begin(); i != vec.end(); i++) {

                //get the global vertex to which the global edge points and find the local vertex holding this
                //global one
                GlobalEdge global = global_extractor()(*i);
                GlobalVertex target = (global.first == fusion::at_c<2>(cg[nv])) ? global.second : global.first;
                std::pair<LocalVertex, bool> res = cg.getContainingVertex(target);
                //if(!res.second) TODO: throw

                //get or create the edge between the new vertex and the target
                LocalEdge e;
                bool done;
                boost::tie(e,done) = boost::edge(nv, res.first, cg);
                if (!done ) boost::tie(e,done) = boost::add_edge(nv, res.first, cg);
                //if(!done) TODO: throw

                //push the global edge to the local edge
                cg[e].push_back(*i);
            };
        }

        //all global edges concerning the move vertex are processed and it is moved to the subcluster,
        //lets destroy it in the local cluster
        boost::clear_vertex(v, *this);
        boost::remove_vertex(v, *this);

        return nv;
    };

    /**
     * @brief Move a vertex to the parent cluster.
     *
     * This function moves a vertex one step up in the subcluster hierarchie and reconnects all other vertices and clusters.
     * The moved vertex will hold it's global descriptor but get a new local one assigned (the one returned). The same
     * stands for all edges which use the moved vertex: global descriptors stay the same, but they are moved to new
     * local edges. Note that this function is the inverse of vertexToSubcluster, and doing Pseudocode:
     * vertexToParent(vertexToSubcluster(v)) does nothing (only the local descriptor of the moved vertex is
     * diffrent afterwards). 
     * 
     * @param v Local vertex which should be moved to the parents cluster
     * @return LocalVertex the local descriptor of the moved vertex, valid in the parent cluster only.
     **/
    LocalVertex vertexToParent( LocalVertex v ) {

        //if(isRoot()) TODO:throw

        //create new vertex
        vertex_bundle& vb = (*this)[v];
        LocalVertex nv = boost::add_vertex(vb, parent());
        GlobalVertex gv= fusion::at_c<2>(vb);
	
        //get all out_edges of this cluster in the parentcluster (because only they can hold relevant global_Edgs)
        std::vector<LocalEdge> edge_vec;
        LocalVertex this_v = parent().getClusterVertex(*this);
        std::pair<local_out_edge_iterator, local_out_edge_iterator> it = boost::out_edges(this_v, parent());
        for (; it.first != it.second; it.first++) {
            //iterate all global edges and find relevant ones
            edge_bundle& vec = parent()[*it.first];
            typename edge_bundle::iterator i = vec.begin();
            while ( i != vec.end() ) {

                GlobalEdge global = global_extractor()(*i);
                GlobalVertex target;
                if (global.first == gv) target = global.second;
                else if (global.second == gv) target = global.first;
                else {
                    i++;
                    continue;
                }

                std::pair<LocalVertex, bool> res = parent().getContainingVertex(target);

                //get or create the edge between the new vertex and the target
                LocalEdge e;
                bool done;
                boost::tie(e,done) = boost::edge(nv, res.first, parent());
                if (!done ) boost::tie(e,done) = boost::add_edge(nv, res.first, parent());
                //if(!done) TODO: throw

                //push the global edge bundle to the new local edge and erase it in the old
                parent()[e].push_back(*i);
                i = vec.erase(i);
            }
            //see if we should destroy this edge (no global edges remain in local one)
            if (vec.empty()) edge_vec.push_back(*it.first);
        }

        //create a edge between new vertex and this cluster and add all global edges from within this cluster
        it = boost::out_edges(v, *this);
        LocalEdge e;
        if (it.first != it.second) e = boost::add_edge(nv, this_v, parent()).first;
        for (; it.first != it.second; it.first++) {
            edge_bundle& ep = (*this)[*it.first];
            edge_bundle& nep = parent()[e];
            nep.insert(nep.end(), ep.begin(), ep.end());
        }

        //all global edges concerning the move vertex are processed and it is moved to the parent,
        //lets destroy it in the local cluster
        boost::clear_vertex(v, *this);
        boost::remove_vertex(v, *this);

        //it's possible that some local edges in the parent are empty now, let's destroy them
        for (std::vector<LocalEdge>::iterator it=edge_vec.begin(); it!=edge_vec.end(); it++)
            boost::remove_edge(*it, parent());
	
	return nv;
    };


    /********************************************************
    * Stuff
    * *******************************************************/

protected:
    ClusterGraph* m_parent;
    ClusterMap	  m_clusters;
    details::IDpointer 	  m_id;

    std::pair<LocalVertex, bool> getContainingVertex(GlobalVertex id, bool recursive = true) {

        //check all vertices if they are the id
        std::pair<local_vertex_iterator, local_vertex_iterator>  it = boost::vertices(*this);
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

    template<typename functor>
    typename functor::result_type apply_to_bundle(LocalVertex k, functor f) {
        return f((*this)[k]);
    };
    template<typename functor>
    typename functor::result_type apply_to_bundle(LocalEdge k, functor f) {
        return f((*this)[k]);
    };

    template<typename functor>
    typename functor::result_type apply_to_bundle(GlobalVertex k, functor f) {

        //check all vertices if they are the id
        std::pair<local_vertex_iterator, local_vertex_iterator>  it = boost::vertices(*this);
        for (;it.first != it.second; it.first++) {
            vertex_prop& p = (*this)[*it.first];
            if ( k == fusion::at_c<2>(p) )
                return f(p);
        }

        //check all clusters if they have the object
        for (cluster_iterator it = m_clusters.begin(); it != m_clusters.end(); it++) {
            return (*it.second)->apply_to_bundle<functor>(k, f);
        }

        // return typename functor::base_type(); TODO: Throw (propeties return reference, but cant init a reference temporarily)
    };

    template<typename functor>
    typename functor::result_type apply_to_bundle(GlobalEdge k, functor f) {

        LocalVertex v1,v2;
        bool d1,d2;
        boost::tie(v1,d1) = getContainingVertex(k.first);
        boost::tie(v2,d2) = getContainingVertex(k.second);

        if ( d1&&d2 ) {
            if ( (v1==v2) && isCluster(v1) ) return m_clusters[v1]->apply_to_bundle(k, f);
            else {
                LocalEdge e;
                bool done;
                boost::tie(e, done) = boost::edge(v1,v2,*this);
                //if(!done) TODO: throw, as there has to be a edge!
                return f((*this)[e]);
            };
        }
        //return std::make_pair(NULL, NULL); TODO:Throw
    };
};

} //namespace solver

#endif // CLUSTERGRAPH_HPP
