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

#ifndef CLUSTERGRAPH_HPP
#define CLUSTERGRAPH_HPP

#include <map>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/copy.hpp>
#include <boost/graph/iteration_macros.hpp>

#include <boost/mpl/transform.hpp>
#include <boost/mpl/find.hpp>

#include <boost/iterator/transform_iterator.hpp>
#include <boost/iterator/filter_iterator.hpp>

#include <boost/fusion/include/as_vector.hpp>
#include <boost/fusion/include/mpl.hpp>
#include <boost/fusion/container/vector.hpp>

#include <boost/fusion/include/algorithm.hpp>
#include <boost/fusion/include/at.hpp>
#include <boost/fusion/include/at_c.hpp>
#include <boost/fusion/include/make_vector.hpp>

#include <boost/bind.hpp>

#include "opendcm/core/property.hpp"

#include <Eigen/Core>

namespace mpl = boost::mpl;
namespace fusion = boost::fusion;

namespace dcm {
namespace graph {
    
using namespace details;

/** @addtogroup Core
 * @{
 * */

/** @addtogroup ClusterGraph
 * @{*/

//we need a way to store a pointer to a graph in a type independend way
struct ClusterGraphBase {};

//type of a graph cluster
/**
 * @brief Add a type to clusters
 *
 * Allows to specify special types to ClusterGraphs and make a it possibe to distuingish between
 * diffrent purposes. The cluster types need to be int.
 **/
struct type_info {
    //states the type of a cluster
    typedef int type;
};
//cluster in graph changed?
/**
 * @brief Was the cluster changed?
 *
 * Adds a boolean to the cluster which indicates if the cluster was changed since the last
 * processing. It should be set to true if vertices and edges were added or removed, Subclusters
 * created or deleted and so on.
 **/
struct changed {
    typedef bool type;
};
/**
 * @brief Add an index to vertices
 *
 * Most boost graph algorithms need a index for vertices, ranging from 0 to vertex count. As this can
 * be useful for many other things it is added as vertex property.
 **/
struct vertex_index {
    typedef int type;
};
/**
 * @brief Add an index to edges
 *
 * Most boost graph algorithms need a index for edges, ranging from 0 to edge count. As this can
 * be useful for many other things it is added as edge property.
 **/
struct edge_index {
    typedef int type;
};


typedef boost::adjacency_list_traits<boost::listS, boost::listS, boost::undirectedS> list_traits;


/**
 * @brief A type to be used as identifier for vertices and edges
 *
 * Vertices and edges need to be identified in a stable(safe/load), unique(over multiple clusters) and
 * comparable manner. The bgl vertex and edge discriptors don't fullfill this need as they have a direct
 * relation to the graphs storage. Therefore they change value on moving entitiys to diffrent clusters or
 * clone actions. This class is used to overcome this problem.
 **/
typedef int universalID;

/**
 * @brief Generator for unique identifiers
 *
 * The universalID used to identify vertices and edges globaly need to be unique and therefore can't be
 * created at good will. This generator creates universalID's in a incremental manner and is intended to
 * to be shared between all graphs of a system, so that all created ID's are unique.
 **/
struct IDgen {
    universalID* counter;

    IDgen() {
        counter = new universalID(10);
    };
    IDgen(universalID id) {
        counter = new universalID(id);
    };
    ~IDgen() {
        delete counter;
    };
    /**
     * @brief Generates a new unique ID
     *
     * @return :universalID
     **/
    universalID generate() {
        return ++ (*counter);
    };
    /**
     * @brief Returns the amount if generated ID's
     *
     * As universalID's are integers the returned count is a ID and can therefore also be used as the last
     * created ID.
     *
     * @return :universalID
     **/
    universalID count() {
        return (*counter);
    };
    /**
     * @brief Set the current value for incremental creation
     *
     * ID's are created incrementaly and if a specific startingpoint is wished it can be set here by
     * supplying the last created ID or the amount of totaly created ID's
     *
     * @param id The last created ID
     * @return void
     **/
    void setCount(universalID id) {
        *counter = id;
    };
};

/**
 * @brief Exception thrown from the graph at any occuring error
 **/
struct cluster_error : virtual boost::exception {};

/**
 * @brief Pointer type to share a common ID generator @ref IDgen
 **/
typedef std::shared_ptr<IDgen> IDpointer;


/** @name Descriptors */
/**@{
 * @brief Identifier for local vertices
 *
 * The boost graph library works with identifiers for vertices which directly relate to there storage.
 * Therefore they can be used only in the relevant cluster, they are local. These are the descriptors
 * which need to be used for all bgl algorithms.
 **/
typedef list_traits::vertex_descriptor 	LocalVertex;

/**
 * @brief Indentifier for local edge
 *
 * The boost graph library works with identifiers for edges which directly relate to there storage.
 * Therefore they can be used only in the relevant cluster, they are local. These are the descriptors
 * which need to be used for all bgl algorithms.
 **/
typedef list_traits::edge_descriptor 		LocalEdge;

/**
 * @brief Identifier for global vertex
 *
 * To overcome the locality of the bgl vertex descriptors a global alternative is introduced. This descriptor
 * is unique over clusters and stable on moves and clones.
 **/
typedef universalID 				GlobalVertex;

/**
 * @brief Identifier for global edge
 *
 * To overcome the locality of the bgl edge discriptors a global alternative is introduced. This descriptor
 * is unique over clusters and stable on moves and clones. It holds it's source and target also as global
 * descriptors of type GlobalVertex and has a unique ID in form of a universalID assigned.
 **/
struct 	GlobalEdge {
    GlobalVertex source;
    GlobalVertex target;
    universalID ID;

    bool operator== (const GlobalEdge& second) const {
        return ID == second.ID;
    };
    bool operator!= (const GlobalEdge& second) const {
        return ID != second.ID;
    };
    bool valid() {
        return ID > 9;
    };
};

/**
 * @brief Store a global vertex as proeprty
 *
 * Convienience property to store a global vertex
 **/
struct VertexProperty {
    typedef GlobalVertex type;
    struct  change_tracking{};
};
/**
 * @brief Store a global edge as property
 *
 * Convienience property to store a global edge
 **/
struct EdgeProperty {
    typedef GlobalEdge type;
    struct  change_tracking{};
};

/**
 * @brief Store multiple global edges as property
 *
 * Convienience property to store multiple global edges in a std::vector
 **/
struct MultiEdgeProperty {
    typedef std::vector<GlobalEdge> type;
    struct  change_tracking{}; 
};
/**@}*/

//Define vertex and edge properties which are always added for use in the boost graph library algorithms
//or which are needed in the internal algorithms
typedef mpl::vector2<vertex_index, VertexProperty> bgl_v_props;
typedef mpl::vector1<edge_index> bgl_e_props;
typedef mpl::vector1<EdgeProperty> bgl_ge_props;

//allow to hold multiple global edges in a property
template<typename ge_props>
struct GlobalEdgeProperty {
    typedef typename ensure_properties<ge_props, bgl_ge_props>::type props;
    typedef std::vector< PropertyOwner<props> > type;
    struct change_tracking{};
};


/**
 * @brief A graph that can be stacked in a tree-like manner without loosing it connections
 *
 * This is basicly a boost adjacency_list with single linked lists 'listS' as storage for vertices and
 * edges. The edges are undirected. This allows to use all boost graph algorithms and provides therefore
 * an comprehensive way for analysing and manipulating its content. It further extends the class with the
 * possibility to cluster its content and to add properties to all entitys. For more
 * information, see the module ClusterGraph
 *
 * @tparam edge_prop a mpl::vector with properties which are added to local edges
 * @tparam globaledge_prop a mpl::vector with properties which are added to global edges
 * @tparam vertex_prop a mpl::vector with properties which are added to vertices
 * @tparam cluster_prop a mpl::vector with properties which are added to all clusters
 **/
template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
class ClusterGraph :
    public boost::adjacency_list < boost::listS, boost::listS, boost::undirectedS,
    PropertyOwner< typename ensure_properties<vertex_prop, bgl_v_props>::type >,
    PropertyOwner< typename ensure_properties<edge_prop, 
                                typename mpl::push_back<bgl_e_props,GlobalEdgeProperty<globaledge_prop>>::type>::type>>,
    public PropertyOwner<typename ensure_property<cluster_prop, changed>::type>,
    public ClusterGraphBase,
    public boost::noncopyable,
    public std::enable_shared_from_this<ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop> > {

public:
    
    /**
     * @brief Convinience type for the GlobalEdgeProperty
     * 
     * Just a typedef to make using the GlobalEdgeProperty od the local edge simpler
     */
    typedef GlobalEdgeProperty<globaledge_prop> GEdgeProperty;
    
    /**
     * @brief mpl::vector with all edge properties
     *
     * The edge properties supplied as template argument to the ClusterGraph are extended with graph
     * specific properties, for example a edge_index_prop. These extra properties are intendet to be
     * used with boost graph algorithms as property maps. They need to be in specefied by the ClusterGraph
     * as they are used within it's implementation. If the graph specific properties are already a part
     * of the given property sequence, nothing happens, they are not added twice.
     **/
    typedef typename ensure_properties<edge_prop, typename mpl::push_back<bgl_e_props,
        GlobalEdgeProperty<globaledge_prop>>::type>::type   edge_properties;

    /**
     * @brief mpl::vector with all global edge properties
     *
     * The global edge properties supplied as template argument to the ClusterGraph are only exposed to ensure
     * a consistent interface. Global edges are not used by any graph algorithm, so we don't need any properties
     * here
     **/
    typedef typename ensure_properties<globaledge_prop, bgl_ge_props>::type   globaledge_properties;

    /**
     * @brief mpl::vector with all vertex properties
     *
     * The vertex properties supplied as template argument to the ClusterGraph are extended with graph
     * specific properties as vertex_index_prop. These extra properties are intendet to be
     * used with boost graph algorithms as property maps. They need to be in specefied by the ClusterGraph
     * as they are used within it's implementation.If the graph specific properties are already a part
     * of the given property sequence, nothing happens, they are not added twice.
     **/
    typedef typename ensure_properties<vertex_prop, bgl_v_props>::type vertex_properties;

    /**
     * @brief The property bundle for GlobalEdges
     *
     * A local edge in a cluster can hold multiple gloabal ones. Therefor we need an extra bundle for
     * the GlobalEdges. This bundle holds the objects which are added to that global edge and it's identifier.
     * Note that global edges don't have properties, these are only for local ones.
     **/
    typedef PropertyOwner<globaledge_properties> edge_bundle_single;
    /**
     * @brief The property bundle for local edges
     *
     * Local edges can hold multiple global ones, we therefore need a std::vector of global edges. As
     * they are fully described by a edge_bundle_single we store those. Also local edges can have properties,
     * so store a fusion sequence of them too.
     **/
    typedef PropertyOwner<edge_properties> edge_bundle;
    /**
     * @brief Iteator to access all edge_bundle_single stored in a edge_bundle
     **/
    typedef typename std::vector< edge_bundle_single >::iterator edge_single_iterator;
    /**
     * @brief Property bundle for local vertices
     *
     * This bundle is simpler than the edge one, as every vertex has on single bundle. We therefore
     * store the global descriptor for identification andthe fusion sequence with the properties all in one
     * bundle.
     **/
    typedef PropertyOwner<vertex_properties> vertex_bundle;

    /**
     * @brief The adjacency_list type ClusterGraph inherited from
     **/
    typedef boost::adjacency_list < boost::listS, boost::listS,
            boost::undirectedS, vertex_bundle, edge_bundle > Graph;

    typedef std::enable_shared_from_this<ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop> > sp_base;

    //if changed_prop is not a property we have to add it now
    typedef typename ensure_property<cluster_prop, changed>::type cluster_properties;

    typedef typename boost::graph_traits<Graph>::vertex_iterator   local_vertex_iterator;
    typedef typename boost::graph_traits<Graph>::edge_iterator     local_edge_iterator;
    typedef typename boost::graph_traits<Graph>::out_edge_iterator local_out_edge_iterator;

    typedef std::map<LocalVertex, std::shared_ptr<ClusterGraph> > ClusterMap;


    struct global_extractor  {
        typedef GlobalEdge result_type;
        template<typename T>
        result_type operator()(T& bundle) const;
    };

    struct global_vertex_extractor  {
        typedef GlobalVertex result_type;
        ClusterGraph& graph;
        global_vertex_extractor(ClusterGraph& g);
        result_type operator()(LocalVertex& v) const;
    };

    //iterators
    /**
     * @brief Iterator for global edge descriptors \ref GlobalEdge
     **/
    typedef boost::transform_iterator<global_extractor, edge_single_iterator> global_edge_iterator;

    /**
     * @brief Iterator for global vertex descriptor \ref GlobalVertex
     **/
    typedef boost::transform_iterator<global_vertex_extractor, local_vertex_iterator> global_vertex_iterator;

    /**
     * @brief Iterator for clusters
     *
     * Allows to iterate over all subclusters.
     **/
    typedef typename ClusterMap::iterator cluster_iterator;
    /**
     * @brief Const equivalent to \ref cluster_iterator
     **/
    typedef typename ClusterMap::const_iterator const_cluster_iterator;

    /**
     * @brief Basic constructor
     *
     * This constructor creates a empty cluster with a new ID generator. This is to be used on initial
     * clustergraph creation, so only for the very first cluster.
     **/
    ClusterGraph() : m_id(new IDgen) {};

    /**
     * @brief Dependent constructor
     *
     * This constructor creates a new cluster, but uses the given cluster as parent. It will therefore
     * create a tree-like relationship. Be aware, that the new cluster is not added to the parents
     * subcluster list, that has to be done manualy. The new cluster shares the parents ID generator.
     *
     * @param g the parent cluster graph
     **/
    ClusterGraph(std::shared_ptr<ClusterGraph> g) : m_parent(g), m_id(new IDgen) {
        if(g)
            m_id = g->m_id;
    };

    ~ClusterGraph() {};

    /**
     * @brief Copys the Clustergraph into a new one
     *
     * Copys this cluster and all subclusters into the give one, which is cleared bevore copying. Be
     * aware that all properties are only copied, and as some may hold shared pointers you may have to
     * clone them. If needed this can be done with the supplied functor, which receives all copied objects
     * to his function operator which returns the new value.
     * @param into The Graph that should be a copy of this
     * @param functor The function objects which gets the graph properties and returns the ones for the
     * copied graph
     */
    template<typename Functor>
    void copyInto(std::shared_ptr<ClusterGraph> into, Functor& functor) const;

    /**
     * @brief Compare by adress, not by content
     * @param other the cluster to compare with
     * @return bool if this is the same cluster in memory
     **/
    template<typename T>
    bool operator== (const T& other) const;

    /**
     * @brief Compare by adress, not by content
     * @param other the cluster to compare with
     * @return bool if this is the not same cluster in memory
     **/
    template<typename T>
    bool operator!= (const T& other) const;

    /**
     * @brief Set diffrent behaviour for changed markers
     *
     * Some methods of the ClusterGraph set it's changed_prop to true. Thats sensible, as they change
     * the graph. However, there are situations where you want to use the methods but don't want the change
     * marked. For example recreations while cloning. This method can be used to disable the changed setting.
     * @param on Turn change markers on or of
     * @return void
     **/
    void setCopyMode(bool on);

    //Make sure the compiler finds the base class setters even with equal named functions in this class
    using PropertyOwner<cluster_properties>::getProperty;
    using PropertyOwner<cluster_properties>::setProperty;
    using PropertyOwner<cluster_properties>::isPropertyChanged;
    using PropertyOwner<cluster_properties>::acknowledgePropertyChange;
    using PropertyOwner<cluster_properties>::hasPropertyChanges;
    using PropertyOwner<cluster_properties>::acknowledgePropertyChanges;

    /**
     * @brief Set the property of a owned cluster
     *
     * Makes it easy to set a property of a subcluster without retrieving it first
     *
     * @tparam P the property type which shall be set
     * @param v the local vertex which describes the subcluster
     **/
    template<typename P>
    typename P::type& getSubclusterProperty(LocalVertex v);

     /**
     * @brief Creates filtered iterator ranges dependend on predicate
     * 
     * Often you want only iterate over a subset of graph elements dependend on some kind of predicate. This
     * function allows you to specify a user defined predicate and returns new iterators which represent only 
     * the elements of the input range for which the predicate returns true. You can easily construct such
     * filter iterators yourself, but this function reduces the amount of code needed. The predicate needs to 
     * be a function object which gets the graph object as constructor parameter. The operator() is called with
     * the objects the iterator type gives when dereferenced.
     * 
     * Don't use this method when you want to filter by properties and then access those properties, as this 
     * means to access them twice and can be done more efficient.
     * 
     * @tparam Predicate a predicate function object which returns true or false, dependend if the given
     * 		object shall be used in the sequence
     * @param range a std::pair of iterators as returned by boos graph iteration functions
     * @return std::pair of new iterators which represent the a subset of range
     */
    template<typename Predicate, typename Iterator>
    std::pair<boost::filter_iterator<Predicate, Iterator>, boost::filter_iterator<Predicate, Iterator> >
    filterRange(std::pair<Iterator, Iterator> range);

    /* *******************************************************
     * Subclustering
     * *******************************************************/

    /**
     * @brief Creates a new subcluster
     *
     * As clusters can be stacked in a tree like manner, this function can be used to create new
     * children. It automaticly adds it to the subcluster list and adds it to the graph. The new
     * subcluster is fully defined by its object and the vertex descriptor which is it's position
     * in the current cluster.
     *
     * @return :pair< std::shared_ptr< ClusterGraph >, LocalVertex > Subcluster and its descriptor
     **/
    std::pair<std::shared_ptr<ClusterGraph>, LocalVertex> createCluster();

    /**
     * @brief Returns the parent cluster
     *
     * In the stacked cluster hirarchy most clusters have a parent whcih can be accessed with this function.
     * However, the toplevel cluster dos nothave a parent and a empty shared_ptr is returned.
     *
     * @return :shared_ptr< ClusterGraph > the parent cluster or empty pointer
     **/
    std::shared_ptr<ClusterGraph> parent();

    /**
     * @brief const version of \ref parent()
     *
     * @return :shared_ptr< ClusterGraph >
     **/
    const std::shared_ptr<ClusterGraph> parent() const;

    /**
     * @brief Is this the toplevel cluster?
     *
     * @return bool if it is
     **/
    bool isRoot() const;

    /**
     * @brief Returns the toplevel cluster
     *
     * @return :shared_ptr< ClusterGraph >
     **/
    std::shared_ptr<ClusterGraph>	 root();

    /**
     * @brief const equivalent of \ref root()
     *
     * @return :shared_ptr< ClusterGraph >
     **/
    const std::shared_ptr<ClusterGraph> root() const;

    /**
     * @brief Iterators for all subclusters
     *
     * A pair with two \ref cluster_iterator is returned which point to the first cluster and
     * to one after the last. #this allows full iteration over all subclusters
     *
     * @return :pair< cluster_iterator, cluster_iterator >
     **/
    std::pair<cluster_iterator, cluster_iterator> clusters();

    /**
     * @brief const equivalent to \ref clusters()
     *
     * @return :pair< const_cluster_iterator, const_cluster_iterator >
     **/
    std::pair<const_cluster_iterator, const_cluster_iterator> clusters() const;

    /**
     * @brief The amount of all subclusters
     *
     * @return :size_t
     **/
    std::size_t numClusters() const;

    /**
     * @brief Check if this vertex is a cluster
     *
     * A subcluster is added as normal vertex to the parent cluster. There is no way to distinguish
     * between clusters and normal vertices with global or local descriptors only. Therefore this
     * function can be used to get information about the type. If it is a cluster, it can be accessed
     * with \ref getVertexCluster
     *
     * @param v The vertex to be checked
     * @return bool is cluster or not
     **/
    bool isCluster(const LocalVertex v) const;

    /**
     * @brief Get the cluster corresponding the discriptor
     *
     * A subcluster is added as normal vertex to the parent cluster. There is no way to access
     * the clusters object with global or local descriptors only. Therefore this
     * function can be used to get the object belonging to the descriptor. If the vertex is not
     * a cluster an empty pointer is returned.
     *
     * @param v The vertex for which the cluster is wanted
     * @return std::shared_ptr<ClusterGraph> the coresponding cluster orempty pointer
     **/
    std::shared_ptr<ClusterGraph> getVertexCluster(LocalVertex v);

    /**
     * @brief Get the vertex descrptor which descripes the clusters position in the graph
     *
     * This function is the inverse to \ref getVertexCluster
     *
     * @param g the graph for which the vertex is searched
     * @return :LocalVertex
     **/
    LocalVertex	getClusterVertex(std::shared_ptr<ClusterGraph> g);

    /**
     * @brief Convinience function for \ref removeCluster
     **/
    template<typename Functor>
    void removeCluster(std::shared_ptr<ClusterGraph> g, Functor& f);
    /**
     * @brief Convinience function for \ref removeCluster
     **/
    void removeCluster(std::shared_ptr<ClusterGraph> g);
    /**
     * @brief Delete all subcluster
     *
     * @return void
     **/

    void clearClusters();

    /**
     * @brief Remove a subcluster and applys the functor to all removed edges and vertices
     *
     * All downstream elements of the local vertex v will be removed after the functor is applied to there
     * edges and vertices. Note that the LocalVertex which represents the cluster to delete is not passed
     * to the functor. When ever the cluster is changed it will be passed to the functor, so that it need
     * to have three overloads: operator()(GlobalEdge), operator()(GlobalVertex), operator()(ClusterGraph&)
     *
     * @param v Local vertex which is a cluster and which should be deleted
     * @param f Functor to apply on all graph elements
     */
    template<typename Functor>
    void removeCluster(LocalVertex v, Functor& f);
    void removeCluster(LocalVertex v);

protected:
    template<typename Functor>
    void remove_vertices(Functor& f, bool recursive = false);


    /* *******************************************************
    * Creation Handling
    * *******************************************************/

public:
    /**
     * @brief Add a vertex to the local cluster
     *
     * @return fusion::vector<LocalVertex, GlobalVertex> the local and global vertex descriptor
     **/
    fusion::vector<LocalVertex, GlobalVertex> addVertex();

    /**
     * @brief Add a vertex to the local cluster with given global identifier
     *
     * Sometimes it is needed to add a vertex with given global identifier. As the global vertex can not
     * be changed after creation, this method can be used to specify the global vertex by which this
     * graph vertex can be identified. The given global vertex is not checked, you need to ensure that
     * it is a unique id or the already existing vertex is returned.
     * The ID generator is changed so that it creates only identifier bigger than v.
     *
     * @return fusion::vector<LocalVertex, GlobalVertex> the local and global vertex descriptor
     **/
    fusion::vector<LocalVertex, GlobalVertex> addVertex(GlobalVertex v);

    /**
     * @brief Iterators of all global vertices in this cluster
     *
     * Returns the iterator for the first global vertex and the end() iterator as reference for
     * iterating
     *
     * @return std::pair< global_vertex_iterator, global_vertex_iterator > global vertex iterators
     **/
    std::pair<global_vertex_iterator, global_vertex_iterator> globalVertices();

    /**
     * @brief Returns the edge between the local vertices
     *
     * This function is the same as boost::edge(source, target, Graph) and only added for convienience.
     *
     * @param source LocalVertex as edge source
     * @param target LocalVertex as edge target
     * @return std::pair<LocalEdge, bool> with the local edge descriptor if existing. The bool value shows if the
     * edge exists or not
     **/
    std::pair<LocalEdge, bool> edge(LocalVertex source, LocalVertex target);
    
    /**
     * @brief Returns the edge between the global vertices
     * 
     * The edge retrievel function for global vertices. Note that the function fails if the global edges are 
     * not part of the same cluster.
     * @param r1 GlobalVertex from which the edge starts
     * @param r2 GlobalVertex where the edge ends
     * @return std::pair<LocalEdge, bool> with the local edge descriptor if existing and a bool value if the edge
     * is valid. It is false if the edge does not exist.
     */
    std::pair<LocalEdge, bool> edge(GlobalVertex source, GlobalVertex target);

    /**
     * @brief Add a edge between two vertices, defined by local descriptors.
     *
     * Add an edge that connects the two vertices and in the local clustergraph and assign the GlobalEdge to it. The
     * LocalVertex parameters should not represent a cluster which would result in the functions failure. If there's
     * already a local edge between the vertices a new global edge will be added and returned. Failure will be
     * recocnisable by a false value in the returned type sequence.
     *
     * @param source The first vertex the edge should connect
     * @param target The second vertex the edge should connect
     * @return fusion::vector<LocalEdge, GlobalEdge, success> with the local and global descriptors of the edge and an bool
     * value indicationg the successful creation.
     **/
    fusion::vector<LocalEdge, GlobalEdge, bool> addEdge(LocalVertex source, LocalVertex target);

    /**
     * @brief Add a edge between two vertices, defined by global descriptors.
     *
     * Adds an edge between vertices which are not nesseccarily in this local cluster and have therefore to be
     * identified with global descriptors. The only condition for source and target vertex is that both must be
     * in the local cluster or any of its subclusters. If thats not the case, the function will fail. On success
     * a new GlobalEdge will be created, but not neccessarily a local one. If the vertices are in different cluster
     * which are already connected the global edge will be added to this connecting local edge. Thats the one returned
     * in the seqence. Note that it's possible that the local edge belongs to another subcluster and therefore can't be
     * used in the local cluster. This case is indicated by the scope return value.
     *
     * @param source The first vertex the edge should connect
     * @param target The second vertex the edge should connect
     * @return fusion:vector< LocalEdge, GlobalEdge, success, scope > with the new global edge descriptor and the local
     * one where it was added. Success indicates if the function was successful and scope shows the validy of the local
     * descriptor in this cluster (true means the edge is in this cluster).
     **/
    fusion::vector<LocalEdge, GlobalEdge, bool, bool> addEdge(GlobalVertex source, GlobalVertex target);

    fusion::vector<LocalEdge, GlobalEdge, bool, bool> addEdgeGlobal(GlobalVertex source, GlobalVertex target);

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
    std::pair<global_edge_iterator, global_edge_iterator> getGlobalEdges(LocalEdge e);

    /**
     * @brief Get the count of all global edges
     *
     * Local edges can hold multiple global ones, for example when they connect at least one cluster. To get the
     * number of all global edges in this local one you can use this function.
     *
     * @param e the local edge for which the global descriptors are wanted
     * @return std::pair<begin, end> with the global_edge_iterator's pointing to the vector<GlobalEdge>'s start
     * and end
     **/
    int getGlobalEdgeCount(LocalEdge e);

    /**
     * @brief Get the local edge which holds the specified global edge.
     *
     * Note that GlobalEdge must be in a local edge of this cluster, means the connected vertices must be in this
     * or one of it's subclusters (but not the same). Also if the containing LocalEdge is not in this cluster, but in one of it's
     * subclusters, the function fails and the returned edge is invalid.
     *
     * @param e GlobalEdge for which the containing local one is wanted
     * @return std:pair< LocalEdge, bool > with the containing LocalEdge and a bool indicator if function was successful.
     **/
    std::pair<LocalEdge, bool> getLocalEdge(GlobalEdge e);

    /**
     * @brief Get the local edge which holds the specified global one and the subcluster in which it is valid.
     *
     * The function only fails when the global edge is hold by a local one upstream in the cluster
     * herarchy.
     *
     * @param e GlobalEdge for which the containing local one is wanted
     * @return fusion::vector<LocalEdge, ClusterGraph*, bool> with the containing LocalEdge, the cluster which holds it and a bool indicator if function was successful.
     **/
    fusion::vector<LocalEdge, ClusterGraph*, bool> getLocalEdgeGraph(GlobalEdge e);

    /**
     * @brief Get the GlobalVertex assiociated with this local one.
     *
     * @param v LocalVertex
     * @return GlobalVertex
     **/
    GlobalVertex getGlobalVertex(LocalVertex v) const;

    /**
     * @brief Get the LocalVertex which corresponds to the golab one
     *
     * The GlobalVertex has to be in this cluster or any of it's subclusters. If its in a subcluster, the returned
     * LocalVertex will represent this cluster. If the GlobalVertex is not in this clusters scope the function fails.
     *
     * @param vertex GlobalVertex for which the local one shall be returned
     * @return std::pair< LocalVertex, bool > The LocalVertex containing the global one and an success indicator
     **/
    std::pair<LocalVertex, bool> getLocalVertex(GlobalVertex vertex);

    /**
     * @brief Get the local vertex which holds the specified global one and the subcluster in which it is valid.
     *
     * The function only fails when the global vertex is hold by a local one upstream in the cluster
     * herarchy.
     *
     * @param v GlobalVertex for which the containing local one is wanted
     * @return fusion::vector<LocalVertex, ClusterGraph*, bool> with the containing LocalVertex, the cluster which holds it and a bool indicator if function was successful.
     **/
    fusion::vector<LocalVertex, std::shared_ptr<ClusterGraph>, bool> getLocalVertexGraph(GlobalVertex v);


    /* *******************************************************
     * Remove Handling
     * *******************************************************/
private:

    template<typename Functor>
    void downstreamRemoveVertex(GlobalVertex v, Functor& f);

    void simpleRemoveEdge(LocalEdge e);


public:
    /**
    * @brief Removes a vertex from the local cluster and applys functor to removed edges
    *
    * Removes the vertex from the local graph and invalidates the global vertex id. Also all edges connecting
    * to this vertex will be removed after the functor was applied to them. The functor needs to implement
    * operato()(GlobalEdge e). Remark that there is no checking done if the vertex is a cluster, so you
    * need to make sure it's not, as removing a clustervertex will not delete the coresponding cluster.
    *
    * @param id Local Vertex which should be removed from the graph
    * @param f functor whose operator(GlobalEdge) is called for every removed edge
    **/
    template<typename Functor>
    void removeVertex(LocalVertex id, Functor& f);
    //no default template arguments for template functions allowed before c++0x, so a little workaround
    void removeVertex(LocalVertex id) ;

    /**
    * @brief Removes a vertex from the cluster or it's subclusters and applys functor to removed edges
    *
    * Removes the vertex from the graph or subclusters and invalidates the global vertex id. Also all edges connecting
    * to this vertex will be removed (upstream and downstream) after the functor was applied to them. The functor
    * needs to implement operato()(LocalEdge edge).
    *
    * @param id Global Vertex which should be removed from the graph
    * @param f functor whose operator(LocalEdge) is called on every removed edge
    **/
    template<typename Functor>
    void removeVertex(GlobalVertex id, Functor& f);
    //no default template arguments for template functions allowed before c++0x, so a little workaround
    void removeVertex(GlobalVertex id);

    /**
    * @brief Removes a global Edge from the cluster or it's subclusters
    *
    * Removes the edge from the graph or subclusters and invalidates the global edge id. If the local edge holds
    * only this global one it will be removed also.
    *
    * @param id Global Edge which should be removed from the graph
    * @return bool indicates if the global id could be removed
    **/
    void removeEdge(GlobalEdge id);

    /**
    * @brief Removes a local edge from the cluster and calls the functor for all removed global edges
    *
    * Removes the edge from the graph and invalidates the global edges. The Functor needs to provide
    * operator()(GlobalEdge). If no functor is needed just use boost::remove_edge.
    *
    * @param id Global Edge which should be removed from the graph
    * @param f functor whoms operator(GlobalEdge) is called
    * @return bool indicates if the global id could be removed
    **/
    template<typename Functor>
    void removeEdge(LocalEdge id, Functor& f);


    /* *******************************************************
     * Property Handling
     * *******************************************************/

    /**
    * @brief Get the desired property at the specified vertex or edge
    *
    * This function allows to access the properties stored in the graph. If no property of the desired type
    * was set before, a default construced will be returned. Accessing the property at a global edge will return
    * the property of the holding local edge.
    *
    * @remark This function is reentrant and can be called with the same data from multiple threads
    * 
    * @tparam property the property type which shall be returned
    * @param k local or global Vertex/Edge descriptor for which the property is desired
    * @return property::type& the reference to the desired property
    **/
    template<typename property, typename key>
    const typename property::type& getProperty(key k);

    /**
     * @brief Set a property at the specified vertex or edge
     *
     * Sets the given value at the given key. Note that every entity can hold one of each property
     * 
     * @remark The function is reentrant but can not be called on the same data from multiple threads
     *
     * @tparam property the property type which shall be set
     * @param k local or global Vertex/Edge descriptor for which the property should be set
     * @param val the property value which should be stored
     * @return void
     **/
    template<typename property, typename key>
    void setProperty(key k, const typename property::type& val);

    /**
     * @brief Set a property at the specified vertex or edge
     *
     * Sets the given value at the given key. Note that every entity can hold one of each property
     *
     * @tparam property the property type which shall be queried
     * @param k local or global Vertex/Edge descriptor which should be accessed
     * @return bool true if the proeprty was changed after the last acknowledgement
     **/
    template<typename property, typename key>
    bool isPropertyChanged(key k);

    /**
     * @brief Acknowledge property change
     *
     * Marks the property as unchanged. This can be used to notice that the change was processed.
     *
     * @tparam property the property type which shall be queried
     * @param k local or global Vertex/Edge descriptor which should be accessed
     * @return bool true if the proeprty was changed after the last acknowledgement
     **/
    template<typename property, typename key>
    void acknowledgePropertyChange(key k);

    /**
     * @brief Check if any property was changed
     *
     * Checks if any property of was changed. Note that this function only checks the local edge if a
     * local edge descriptor is given, not the global ones. To check for changes in a local edge and all its 
     * global edges use \ref edgeChanged 
     * 
     * @param k local or global Vertex/Edge descriptor which should be accessed
     * @return bool true if any property has the change flag set to true, i.e. isPropertyChanged() returns true
     */
    template<typename key>
    bool hasPropertyChanges(key k);

    /**
     * @brief Acknowledge every property
     *
     * Sets the change flag for every property to false. Note that this function only acknowledges the local edge
     * if a local edge descriptor is given, not the global ones. To check for changes in a local edge and all its 
     * global edges use \ref acknowledgeEdgeChanges 
     * @param k local or global Vertex/Edge descriptor which should be accessed
     */
    template<typename key>
    void acknowledgePropertyChanges(key k);
    
    /**
     * @brief Check if the local edge has seen any changes
     * 
     * This function checks the edge for changes, this means it checks for added/removed global edges, changed 
     * properties and changed global edges.
     * 
     * @param e Local edge to check
     * @return bool true if the changes happend 
     */
    bool edgeChanged(LocalEdge e);
    
    /**
     * @brief Set the edge to unchanged
     * 
     * Markes all edge properties and all global edge properties as unchanged and leads to edgeChanged call for
     * this edge to return false.
     * @remark The function is reentrant but is not save to be called with the same edge from different threads
     *  
     * @param e Local edge to acknowledge
     * @return void
     */
    void acknowledgeEdgeChanges(LocalEdge e);

    /**
     * @brief recreate the internal index maps for edges and vertices
     *
     * Quite many boost graph algorithms need the indices for vertices and edges which are provided by property
     * maps. As we use list, and not vector, as underlaying storage we don't get that property for free and
     * need to create it ourself. To ease that procedure the internal property vertex_index_prop and edge_index_prop
     * can be used as property maps and can be initialized by calling this function.
     *
     * @return void
     **/
    void initIndexMaps();

    //possible predicates based on properties to filter the graph iteration
    //can be used together with filter_iterator
    template<typename Property>
    struct property_changed {
        std::shared_ptr<ClusterGraph> graph;
        property_changed(std::shared_ptr<ClusterGraph> g) : graph(g) {};

        template<typename It>
        bool operator()(It it) {
            return graph->template isPropertyChanged<Property>(it);
        };
    };

    template<typename Property>
    struct property_unchanged {
        std::shared_ptr<ClusterGraph> graph;
        property_unchanged(std::shared_ptr<ClusterGraph> g) : graph(g) {};

        template<typename It>
        bool operator()(It it) {
            return !(graph->template isPropertyChanged<Property>(it));
        };
    };

    struct property_changes {
        std::shared_ptr<ClusterGraph> graph;
        property_changes(std::shared_ptr<ClusterGraph> g) : graph(g) {};

        template<typename It>
        bool operator()(It it) {
            return graph->template hasPropertyChanges(it);
        };
    };
    
    template<typename Property, typename Property::type value>
    struct property_value {
        std::shared_ptr<ClusterGraph> graph;
        property_value(std::shared_ptr<ClusterGraph> g) : graph(g) {};

        template<typename It>
        bool operator()(It it) {
            return (graph->template getProperty<Property>(it) == value);
        };
    };
    
    struct edge_changed {
        std::shared_ptr<ClusterGraph> graph;
        edge_changed(std::shared_ptr<ClusterGraph> g) : graph(g) {};

        template<typename It>
        bool operator()(It it) {
            return graph->template edgeChanged(it);
        };
    };
    
    struct edge_unchanged {
        std::shared_ptr<ClusterGraph> graph;
        edge_unchanged(std::shared_ptr<ClusterGraph> g) : graph(g) {};

        template<typename It>
        bool operator()(It it) {
            return !(graph->template edgeChanged(it));
        };
    };


    /********************************************************
    * Vertex and Cluster moving
    * *******************************************************/

    /**
     * @brief Move a vertex to a subcluster
     *
     * Overloaded convinience function which fetches the local descriptor for the cluster reference and calls
     * the full parameter equivalent. Both cluster and vertex must be in the local cluster.
     *
     * @param v the LocalVertex to be moved
     * @param cg reference to the subcluster to which v should be moved
     * @return LocalVertex the local descriptor of the moved vertex in the subcluster
     **/
    LocalVertex moveToSubcluster(LocalVertex v, std::shared_ptr<ClusterGraph> cg);

    /**
     * @brief Move a vertex to a subcluster
     *
     * Overloaded convinience function which fetches the the cluster reference for the local descriptor and calls
     * the full parameter equivalent. Both cluster and vertex must be in the local cluster.
     *
     * @param v the LocalVertex to be moved
     * @param Cluster the local vertex descriptor representing the subcluster to which v should be moved
     * @return LocalVertex the local descriptor of the moved vertex in the subcluster
     **/
    LocalVertex moveToSubcluster(LocalVertex v, LocalVertex Cluster);

    /**
     * @brief Move a vertex to a subcluster
     *
     * This function moves the LocalVertex to the subcluster and reconnects all other vertices and clusters. The
     * moved vertex will hold it's global descriptor but get a new local one assigned (the one returned). The same
     * stands for all edges which use the moved vertex: global descriptors stay the same, but they are moved to new
     * local edges. It's allowed to move cluster vertices with this function.
     * The specified cluster has of course to be a valid and direct subcluster, the move vertex also has to be in the
     * local cluster.
     *
     * @param v the LocalVertex to be moved
     * @param Cluster the local vertex descriptor representing the subcluster to which v should be moved
     * @param cg reference to the subcluster to which v should be moved
     * @return LocalVertex the local descriptor of the moved vertex in the subcluster
     **/
    LocalVertex moveToSubcluster(LocalVertex v, LocalVertex Cluster, std::shared_ptr<ClusterGraph> cg);


    /**
     * @brief Move a vertex to the parent cluster.
     *
     * This function moves a vertex one step up in the subcluster hierarchie and reconnects all other vertices and clusters.
     * The moved vertex will hold it's global descriptor but get a new local one assigned (the one returned). The same
     * stands for all edges which use the moved vertex: global descriptors stay the same, but they are moved to new
     * local edges. Note that this function is the inverse of moveToSubcluster, and doing Pseudocode:
     * moveToParent(moveToSubcluster(v)) does nothing (only the local descriptor of the moved vertex is
     * diffrent afterwards).
     *
     * @param v Local vertex which should be moved to the parents cluster
     * @return LocalVertex the local descriptor of the moved vertex, valid in the parent cluster only.
     **/
    LocalVertex moveToParent(LocalVertex v);


    /********************************************************
    * Stuff
    * *******************************************************/

    ClusterMap	  m_clusters;
    int test;

protected:
    std::weak_ptr<ClusterGraph> m_parent;
    IDpointer 	  m_id;
    bool copy_mode; //no changing itself when copying


    /* Searches the global vertex in all local vertices of this graph, and returns the local
     * one which holds the global vertex. If not successfull the local vertex returned will be
     * invalid and the bool parameter will be false. If recursive = true, all subclusters will
     * be seached too, however, if found there the retourned local vertex will be the vertex
     * representing the toplevel cluster holding the global vertex in the initial graph.
     * */
    std::pair<LocalVertex, bool> getContainingVertex(GlobalVertex id, bool recursive = true);

    /* Searches the local vertex holding the specified global one in this and all it's subclusters.
     * If found, the holding local vertex and the graph in which it is valid will be returned.
     * */
    fusion::vector<LocalVertex, std::shared_ptr<ClusterGraph>, bool> getContainingVertexGraph(GlobalVertex id);

    /* Searches the global edge in all local edges of this graph, and returns the local
     * one which holds the global edge. If not successfull the local edge returned will be
     * invalid and the bool parameter will be false.
     * */
    std::pair<LocalEdge, bool> getContainingEdge(GlobalEdge id);

    /* Searches the local edge holding the specified global one in this and all it's subclusters.
     * If found, the holding local edge and the graph in which it is valid will be returned.
     * */
    fusion::vector<LocalEdge, ClusterGraph*, bool> getContainingEdgeGraph(GlobalEdge id);

    template<typename functor>
    typename functor::result_type apply_to_bundle(LocalVertex k, functor f);

    template<typename functor>
    typename functor::result_type apply_to_bundle(LocalEdge k, functor f);

    template<typename functor>
    typename functor::result_type apply_to_bundle(GlobalVertex k, functor f);

    template<typename functor>
    typename functor::result_type apply_to_bundle(GlobalEdge k, functor f);

public:
    //may hold properties which have Eigen3 objects and therefore need alignment
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

/** @} */
/** @} */


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
            res = (edge == e.template getProperty<EdgeProperty>());
        else
            res = (vert == e.template getProperty<EdgeProperty>().source) || (vert == e.template getProperty<EdgeProperty>().target);

        //we are a hit, invoke the functor.
        if(res || vert < 0)
            func(e.template getProperty<EdgeProperty>());

        return res || vert < 0;
    }
};

template<typename prop, typename Graph>
struct get_property_helper {

    typedef typename prop::type base_type;
    typedef const base_type& result_type;
 
    template<typename bundle>
    result_type operator()(bundle& p) {
        return p.template getProperty<prop>();
    }
};

template<typename prop, typename Graph>
struct set_property_helper {

    typedef typename prop::type base_type;
    typedef const base_type& value_type;
    typedef void result_type;
 
    set_property_helper(value_type v) :  value(v) {};

    template<typename bundle>
    result_type operator()(bundle& p) {
        p.template setProperty<prop>(value);
    };

    value_type value;
};

template<typename prop, typename Graph>
struct change_property_helper {

    typedef bool result_type;
 
    template<typename bundle>
    result_type operator()(bundle& p) {
        return p.template isPropertyChanged<prop>();
    }
};

template<typename prop, typename Graph>
struct acknowledge_property_helper {

    typedef void result_type;
 
    template<typename bundle>
    result_type operator()(bundle& p) {
        p.template acknowledgePropertyChange<prop>();
    }
};

template<typename Graph, typename Key>
struct changes_property_helper {

    typedef bool result_type;
 
    template<typename bundle>
    result_type operator()(bundle& p) {
        p.template hasPropertyChanges();
    }
};

template<typename Graph, typename Key>
struct ack_changes_property_helper {

    typedef void result_type;
 
    template<typename bundle>
    result_type operator()(bundle& p) {
        p.template acknowledgePropertyChanges();
    }
};

//Function implementation
//***********************

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
template<typename T>
typename ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::global_extractor::result_type
ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::global_extractor::operator()(T& bundle) const {
    return bundle.getProperty<EdgeProperty>();
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
void ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::copyInto(std::shared_ptr<ClusterGraph> into, Functor& functor) const {

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
        std::shared_ptr<ClusterGraph> ng = std::shared_ptr<ClusterGraph> (new ClusterGraph(into));

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
std::pair<std::shared_ptr< ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop> >, LocalVertex> ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::createCluster() {
    vertex_bundle vp;
    vp.template setProperty<VertexProperty>(m_id->generate());
    LocalVertex v = boost::add_vertex(vp, *this);
    return std::pair<std::shared_ptr<ClusterGraph>, LocalVertex> (m_clusters[v] = std::shared_ptr<ClusterGraph> (new ClusterGraph(sp_base::shared_from_this())), v);
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
inline std::shared_ptr< ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop> >
ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>:: parent()  {
    return std::shared_ptr<ClusterGraph> (m_parent);
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
inline const std::shared_ptr< ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop> >
ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::parent() const     {
    return std::shared_ptr<ClusterGraph> (m_parent);
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
bool ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::isRoot() const {
    return m_parent.expired();
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
std::shared_ptr< ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop> >
ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::root()             {
    return isRoot() ? sp_base::shared_from_this() : parent()->root();
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
const std::shared_ptr< ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop> >
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
std::shared_ptr< ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop> >
ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::getVertexCluster(LocalVertex v) {
    if(isCluster(v))
        return m_clusters[v];

    //TODO:throw if not a cluster
    return std::shared_ptr< ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop> >();//sp_base::shared_from_this();
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
LocalVertex     ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::getClusterVertex(std::shared_ptr<ClusterGraph> g) {
    std::pair<cluster_iterator, cluster_iterator> it = clusters();

    for(; it.first != it.second; it.first++) {
        if((*it.first).second == g)
            return (*it.first).first;
    }

    throw graph::cluster_error() <<  boost::errinfo_errno(12) << error_message("Cluster is not part of this graph");
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
template<typename Functor>
void ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::removeCluster(std::shared_ptr<ClusterGraph> g, Functor& f) {
    removeCluster(getClusterVertex(g), f);
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
void ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::removeCluster(std::shared_ptr<ClusterGraph> g) {
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
        throw graph::cluster_error() <<  boost::errinfo_errno(11) << error_message("Cluster is not part of this graph");

    std::pair<LocalVertex, std::shared_ptr<ClusterGraph> > res = *it;

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
    vp.template setProperty<VertexProperty>(m_id->generate());
    LocalVertex v = boost::add_vertex(vp, *this);

    return fusion::make_vector(v, m_id->count());
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
fusion::vector<LocalVertex, GlobalVertex>
ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::addVertex(GlobalVertex gv) {

    std::pair<LocalVertex, bool> res = getLocalVertex(gv);

    if(!res.second) {
        vertex_bundle vp;
        vp.template setProperty<VertexProperty>(gv);
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
std::pair<LocalEdge, bool>
ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::edge(GlobalVertex source, GlobalVertex target) {
    
    std::pair<LocalVertex, bool> r1 = getLocalVertex(source);
    std::pair<LocalVertex, bool> r2 = getLocalVertex(target);
      
    if(!r1.second || !r2.second)
        return std::make_pair(LocalEdge(), false);
        
    return edge(r1.first, r2.first);
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
    GlobalEdge global = { ((*this) [source]).template getProperty<VertexProperty>(),
                          ((*this) [target]).template getProperty<VertexProperty>(), m_id->generate() };
    edge_bundle_single s;
    s.setProperty<EdgeProperty>(global);
    auto& vec = ((*this) [e]).template getPropertyAccessible<GEdgeProperty>();
    vec.push_back(s);
    ((*this) [e]).template markPropertyChanged<GEdgeProperty>();

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
    s.setProperty<EdgeProperty>(global);
    auto& vec = ((*this) [e]).template getPropertyAccessible<GEdgeProperty>();
    vec.push_back(s);
    ((*this) [e]).template markPropertyChanged<GEdgeProperty>();
    
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

    //note that we can the global edges via the accessible property function as we are sure we do not change it.
    //changes to the global edges itself are tracked by them.
    auto& vec = ((*this) [e]).template getPropertyAccessible<GEdgeProperty>();
    
    global_edge_iterator begin = boost::make_transform_iterator(vec.begin(), global_extractor());
    global_edge_iterator end   = boost::make_transform_iterator(vec.end(), global_extractor());

    return std::pair<global_edge_iterator, global_edge_iterator> (begin, end);
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
int ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::getGlobalEdgeCount(LocalEdge e) {

    return ((*this) [e]).template getProperty<GEdgeProperty>().size();
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
    return ((*this) [v]).template getProperty<VertexProperty>();
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
std::pair<LocalVertex, bool>
ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::getLocalVertex(GlobalVertex vertex) {
    return getContainingVertex(vertex);
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
fusion::vector<LocalVertex, std::shared_ptr< ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop> >, bool>
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
        auto& vec = ((*this) [* (it.first)]).template getPropertyAccessible<GEdgeProperty>();
        vec.erase(std::remove_if(vec.begin(), vec.end(), apply_remove_prediacte<Functor, ClusterGraph> (f, v)), vec.end());
        ((*this) [* (it.first)]).template markPropertyChanged<GEdgeProperty>();
        
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
    auto& vec = ((*fusion::at_c<1> (res)) [fusion::at_c<0> (res)]).template getPropertyAccessible<GEdgeProperty>();
    vec.erase(std::remove_if(vec.begin(), vec.end(), apply_remove_prediacte<placehoder, ClusterGraph> (p, id)), vec.end());
    ((*fusion::at_c<1> (res)) [fusion::at_c<0> (res)]).template markPropertyChanged<GEdgeProperty>();
    
    if(vec.empty())
        boost::remove_edge(fusion::at_c<0> (res), *fusion::at_c<1> (res));    
    
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
template<typename Functor>
void ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::removeEdge(LocalEdge id, Functor& f) {

    auto& vec = ((*this) [id]).template getPropertyAccessible<GEdgeProperty>();
    std::for_each(vec.begin(), vec.end(), boost::bind<void> (boost::ref(apply_remove_prediacte<placehoder, ClusterGraph> (f, -1)), _1));
    ((*this) [id]).template markPropertyChanged<GEdgeProperty>();
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
bool ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::edgeChanged(LocalEdge e) {
  
    if(((*this)[e]).hasPropertyChanges())
        return true;
    
    auto& vec = ((*this) [e]).template getProperty<GEdgeProperty>();
    for(const edge_bundle_single& global : vec) {
        if(global.template hasPropertyChanges())
            return true;
    }
    return false;
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
void ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::acknowledgeEdgeChanges(LocalEdge e) {
  
    ((*this)[e]).acknowledgePropertyChanges();
        
    auto& vec = ((*this) [e]).template getPropertyAccessible<GEdgeProperty>();
    for(edge_bundle_single& global : vec)
        global.acknowledgePropertyChanges();
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
ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::moveToSubcluster(LocalVertex v, std::shared_ptr<ClusterGraph> cg) {

    LocalVertex cv = getClusterVertex(cg);
    return moveToSubcluster(v, cv, cg);
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
LocalVertex
ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::moveToSubcluster(LocalVertex v, LocalVertex Cluster) {

    std::shared_ptr<ClusterGraph> cg = getVertexCluster(Cluster);
    return moveToSubcluster(v, Cluster, cg);
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
LocalVertex
ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::moveToSubcluster(LocalVertex v, LocalVertex Cluster, std::shared_ptr<ClusterGraph> cg) {

    std::pair<local_out_edge_iterator, local_out_edge_iterator> it =  boost::out_edges(v, *this);

    /* add the later removed edges to the coressponding existing edges
     * (or create new edges between adjacent vertices of moved vertex and cluster).
     * also get the edge between cluster and vertex while iterating */
    for(; it.first != it.second; it.first++) {

        LocalVertex target = boost::target(*it.first, *this);
        if(target == v )
            target = boost::source(*it.first, *this);

        if(target != Cluster) {

            //get or create the edge between the old edge target and the cluster
            LocalEdge e;
            bool done;
            boost::tie(e, done) = boost::edge(target, Cluster, *this);

            if(!done)
                boost::tie(e, done) = boost::add_edge(target, Cluster, *this);

            //if(!done) TODO: throw

            auto& ep = ((*this) [*it.first]).template getPropertyAccessible<GEdgeProperty>();
            auto& nep = ((*this) [e]).template getPropertyAccessible<GEdgeProperty>();
            nep.insert(nep.end(), ep.begin(), ep.end());
            ((*this) [e]).template markPropertyChanged<GEdgeProperty>();
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
        auto& vec = ((*this) [moveedge.first]).template getPropertyAccessible<GEdgeProperty>();

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
            auto& gvec = ((*cg) [e]).template getPropertyAccessible<GEdgeProperty>();
            gvec.push_back(*i);
            ((*cg) [e]).template markPropertyChanged<GEdgeProperty>();
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

    GlobalVertex gv = vb.template getProperty<VertexProperty>();

    //get all out_edges of this cluster in the parentcluster (because only they can hold relevant global_Edgs)
    std::vector<LocalEdge> edge_vec;
    LocalVertex this_v = parent()->getClusterVertex(sp_base::shared_from_this());
    std::pair<local_out_edge_iterator, local_out_edge_iterator> it = boost::out_edges(this_v, *parent());

    for(; it.first != it.second; it.first++) {
        //iterate all global edges and find relevant ones
        auto& vec = ((*parent()) [*it.first]).template getPropertyAccessible<GEdgeProperty>();
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
            auto& gvec =  ((*parent()) [e]).template getPropertyAccessible<GEdgeProperty>();
            gvec.push_back(*i);
            ((*parent()) [e]).template markPropertyChanged<GEdgeProperty>();
            
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
        auto& ep  = ((*this) [*it.first]).template getPropertyAccessible<GEdgeProperty>();
        auto& nep = ((*parent()) [e]).template getPropertyAccessible<GEdgeProperty>();
        nep.insert(nep.end(), ep.begin(), ep.end());
        
        ((*parent()) [e]).template markPropertyChanged<GEdgeProperty>();
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
        if(id == ((*this) [*it.first]).template getProperty<VertexProperty>())
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
fusion::vector<LocalVertex, std::shared_ptr< ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop> >, bool>
ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::getContainingVertexGraph(GlobalVertex id) {

    LocalVertex v;
    bool done;
    boost::tie(v, done) = getContainingVertex(id);

    if(!done)
        return fusion::make_vector(LocalVertex(), std::shared_ptr<ClusterGraph>(), false);

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

        if(k == p.template getProperty<VertexProperty>())
            return f(p);
    }

    //check all clusters if they have the object
    fusion::vector<LocalVertex, std::shared_ptr<ClusterGraph>, bool> res = getContainingVertexGraph(k);

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

        //search the global one in the local edge. Note that we can easily use the accessible version 
        //as we dont change the vector. If the functor changes the individual properties we this will be tracked
        //by them individual
        auto& vec = ((*this) [e]).template getPropertyAccessible<GEdgeProperty>();
        edge_single_iterator it;

        for(it = vec.begin(); it != vec.end(); it++) {
            if((*it).template getProperty<EdgeProperty>() == k)
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

} //namespace graph
} //namespace dcm


#endif // CLUSTERGRAPH_HPP





