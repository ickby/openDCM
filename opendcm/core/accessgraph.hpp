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

#ifndef ACCESSGRAPH_HPP
#define ACCESSGRAPH_HPP

#include <map>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/filtered_graph.hpp>
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

/** @addtogroup AccessGraph
 * @{*/

//we need a way to store a pointer to a graph in a type independend way
struct AccessGraphBase {};

//type of a graph cluster
/**
 * @brief Add a type to clusters
 *
 * Allows to specify special types to AccessGraphs and make a it possibe to distuingish between
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
 * @brief Add an index to vertices and edges
 *
 * Most boost graph algorithms need a index for verticesand edges, ranging from 0 to count. As this can
 * be useful for many other things it is added as property.
 **/
struct Index {
    typedef int type;
};


/**
 * @brief The edge/vertex belongs to a certain group
 *
 * It is possible to filter edges and vertices for the group they belong to either by iterating over all filtered 
 * elemente or by using a filtered graph. This allows to apply boost graph algorithms only to a certain subset of 
 * the graph.
 **/
struct Group {
    typedef int type;
    struct default_value {
        int operator()() {
            return 0;
        };
    };
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
 * @brief Exception thrown from the graph at any occuring error
 **/
struct cluster_error : virtual boost::exception {};


/** @name Descriptors */
/**@{
 * @brief Identifier for local vertices
 *
 * The boost graph library works with identifiers for vertices which directly relate to there storage.
 * Therefore they can be used only in the relevant cluster, they are local. These are the descriptors
 * which need to be used for all bgl algorithms.
 **/
typedef list_traits::vertex_descriptor  LocalVertex;

/**
 * @brief Indentifier for local edge
 *
 * The boost graph library works with identifiers for edges which directly relate to there storage.
 * Therefore they can be used only in the relevant cluster, they are local. These are the descriptors
 * which need to be used for all bgl algorithms.
 **/
typedef list_traits::edge_descriptor            LocalEdge;

/**
 * @brief Identifier for global vertex
 *
 * To overcome the locality of the bgl vertex descriptors a global alternative is introduced. This descriptor
 * is unique over clusters and stable on moves and clones.
 **/
typedef universalID                             GlobalVertex;

/**
 * @brief Identifier for global edge
 *
 * To overcome the locality of the bgl edge discriptors a global alternative is introduced. This descriptor
 * is unique over clusters and stable on moves and clones. It holds it's source and target also as global
 * descriptors of type GlobalVertex and has a unique ID in form of a universalID assigned.
 **/
struct  GlobalEdge {
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
typedef mpl::vector3<Index, Group, VertexProperty> bgl_v_props;
typedef mpl::vector2<Index, Group> bgl_e_props;
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
 * information, see the module AccessGraph
 *
 * @tparam edge_prop a mpl::vector with properties which are added to local edges
 * @tparam globaledge_prop a mpl::vector with properties which are added to global edges
 * @tparam vertex_prop a mpl::vector with properties which are added to vertices
 * @tparam cluster_prop a mpl::vector with properties which are added to all clusters
 **/


template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop, 
            template<class, class, class, class, class> class graph_base>
class  AccessGraph : public AccessGraphBase, public boost::noncopyable, 
            public PropertyOwner<typename ensure_property<cluster_prop, changed>::type>,
            public std::enable_shared_from_this<AccessGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop, graph_base> > {

public:
    
    //mare the template parameters accessible from outside
    typedef edge_prop edgeprop;
    typedef globaledge_prop globaledgeprop;
    typedef vertex_prop vertexprop;
    typedef cluster_prop clusterprop;
    
    
    /**
     * @brief Convinience type for the GlobalEdgeProperty
     * 
     * Just a typedef to make using the GlobalEdgeProperty of the local edge simpler
     */
    typedef GlobalEdgeProperty<globaledge_prop> GEdgeProperty;
    
    /**
     * @brief mpl::vector with all edge properties
     *
     * The edge properties supplied as template argument to the AccessGraph are extended with graph
     * specific properties, for example a Index_prop. These extra properties are intendet to be
     * used with boost graph algorithms as property maps. They need to be in specefied by the AccessGraph
     * as they are used within it's implementation. If the graph specific properties are already a part
     * of the given property sequence, nothing happens, they are not added twice.
     **/
    typedef typename ensure_properties<edge_prop, typename mpl::push_back<bgl_e_props,
        GlobalEdgeProperty<globaledge_prop>>::type>::type   edge_properties;

    /**
     * @brief mpl::vector with all global edge properties
     *
     * The global edge properties supplied as template argument to the AccessGraph are only exposed to ensure
     * a consistent interface. Global edges are not used by any graph algorithm, so we don't need any properties
     * here
     **/
    typedef typename ensure_properties<globaledge_prop, bgl_ge_props>::type   globaledge_properties;

    /**
     * @brief mpl::vector with all vertex properties
     *
     * The vertex properties supplied as template argument to the AccessGraph are extended with graph
     * specific properties as Index_prop. These extra properties are intendet to be
     * used with boost graph algorithms as property maps. They need to be in specefied by the AccessGraph
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
     * @brief The adjacency_list type the AccessGraph uses under the hood
     **/
    typedef graph_base < boost::listS, boost::listS,
            boost::undirectedS, vertex_bundle, edge_bundle > Graph;

    typedef std::enable_shared_from_this<AccessGraph<edge_prop, globaledge_prop, vertex_prop, 
                                                     cluster_prop, graph_base> > sp_base;

    //if changed_prop is not a property we have to add it now
    typedef typename ensure_property<cluster_prop, changed>::type cluster_properties;

    typedef typename boost::graph_traits<Graph>::vertex_iterator   local_vertex_iterator;
    typedef typename boost::graph_traits<Graph>::edge_iterator     local_edge_iterator;
    typedef typename boost::graph_traits<Graph>::out_edge_iterator local_out_edge_iterator;

    struct global_extractor  {
        typedef GlobalEdge result_type;
        template<typename T>
        result_type operator()(T& bundle) const;
    };

    struct global_vertex_extractor  {
        typedef GlobalVertex result_type;
        AccessGraph& graph;
        global_vertex_extractor(AccessGraph& g);
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
     * @brief Basic constructor
     *
     * As the access graph does not store the graph data itself the referenced graph needs to be 
     * supplied via this constructor 
     **/
    AccessGraph(Graph& base) : m_graph(base) {};

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
    
protected:
    vertex_bundle& operator[](const LocalVertex& v) {
        return m_graph[v];
    };
    
    edge_bundle& operator[](const LocalEdge& e) {
        return m_graph[e];
    };
    
    
public:                
    
    //boost convinience functions
    Graph& getDirectAccess() {
        return m_graph;
    };
    
    typename boost::graph_traits<Graph>::edges_size_type edgeCount() const {
        return boost::num_edges(m_graph);
    };

    typename boost::graph_traits<Graph>::vertices_size_type vertexCount() const {
        return boost::num_vertices(m_graph);
    };
    
    typename boost::graph_traits<Graph>::degree_size_type outDegree(LocalVertex v) const {
        return boost::out_degree(v, m_graph);
    };
    
    std::pair<local_edge_iterator, local_edge_iterator> edges() {    
        return boost::edges(m_graph);
    };
    
    std::pair<local_vertex_iterator, local_vertex_iterator> vertices() {
        return boost::vertices(m_graph);
    };

    //Make sure the compiler finds the base class setters even with equal named functions in this class
    using PropertyOwner<cluster_properties>::getProperty;
    using PropertyOwner<cluster_properties>::setProperty;
    using PropertyOwner<cluster_properties>::isPropertyChanged;
    using PropertyOwner<cluster_properties>::acknowledgePropertyChange;
    using PropertyOwner<cluster_properties>::hasPropertyChanges;
    using PropertyOwner<cluster_properties>::acknowledgePropertyChanges;

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
     *          object shall be used in the sequence
     * @param range a std::pair of iterators as returned by boos graph iteration functions
     * @return std::pair of new iterators which represent the a subset of range
     */
    template<typename Predicate, typename Iterator>
    std::pair<boost::filter_iterator<Predicate, Iterator>, boost::filter_iterator<Predicate, Iterator> >
    filterRange(std::pair<Iterator, Iterator> range);

   
    /* *******************************************************
    * Creation Handling
    * *******************************************************/

public:
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
    * @brief Access properties non-const at the specified vertex or edge
    *
    * Don't use this unless abselutely nesseccary. It sets the property to changed, no matter if you realy
    * change it or not 8as opposed to the property owner version). This is needed as it is impossible to detect
    * if the reference was changed outside of the owner. Furthermore you should never ever store a refence to a
    * property, as changes can't be tracked either. This function is only available to comply with boost graph
    * property maps and for properties whiche are to big to effieciently be copyed before and after change.
    * @tparam Prop property type which should be accessed
    * @param k local or global Vertex/Edge descriptor for which the property is desire
    * @return Prop::type& a reference to the properties actual value.
    **/
    template<typename Prop, typename key>
    typename Prop::type& getPropertyAccessible(key k);

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
     * need to create it ourself. To ease that procedure the internal property Index_prop and Index_prop
     * can be used as property maps and can be initialized by calling this function.
     *
     * @return void
     **/
    void initIndexMaps();

    //possible predicates based on properties to filter the graph iteration
    //can be used together with filter_iterator
    template<typename Property>
    struct property_changed {
        std::shared_ptr<AccessGraph> graph;
        property_changed(std::shared_ptr<AccessGraph> g) : graph(g) {};

        template<typename It>
        bool operator()(const It it) const {
            return graph->template isPropertyChanged<Property>(it);
        };
    };

    template<typename Property>
    struct property_unchanged {
        std::shared_ptr<AccessGraph> graph;
        property_unchanged(std::shared_ptr<AccessGraph> g) : graph(g) {};

        template<typename It>
        bool operator()(const It it) const {
            return !(graph->template isPropertyChanged<Property>(it));
        };
    };

    struct property_changes {
        std::shared_ptr<AccessGraph> graph;
        property_changes(std::shared_ptr<AccessGraph> g) : graph(g) {};

        template<typename It>
        bool operator()(const It it) const {
            return graph->template hasPropertyChanges(it);
        };
    };
    
    template<typename Property, typename Property::type value>
    struct property_value {
        std::shared_ptr<AccessGraph> graph;
        property_value(std::shared_ptr<AccessGraph> g) : graph(g) {};

        template<typename It>
        bool operator()(const It it) const {
            return (graph->template getProperty<Property>(it) == value);
        };
    };
    
    struct edge_changed {
        std::shared_ptr<AccessGraph> graph;
        edge_changed(std::shared_ptr<AccessGraph> g) : graph(g) {};

        template<typename It>
        bool operator()(const It it) const {
            return graph->template edgeChanged(it);
        };
    };
    
    struct edge_unchanged {
        std::shared_ptr<AccessGraph> graph;
        edge_unchanged(std::shared_ptr<AccessGraph> g) : graph(g) {};

        template<typename It>
        bool operator()(const It it) const {
            return !(graph->template edgeChanged(it));
        };
    };
    
    template<int group>
    struct in_group {
        std::shared_ptr<AccessGraph> graph;
        in_group(std::shared_ptr<AccessGraph> g) : graph(g) {};

        template<typename It>
        bool operator()(const It it) const {
            return (graph->template getProperty<Group>(it) == group);
        };
    };
    
    template<int group>
    struct outof_group {
        std::shared_ptr<AccessGraph> graph;
        outof_group(std::shared_ptr<AccessGraph> g) : graph(g) {};

        template<typename It>
        bool operator()(const It it) const {
            return (graph->template getProperty<Group>(it) != group);
        };
    };


    /********************************************************
    * Stuff
    * *******************************************************/

private:
    Graph& m_graph;

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

//********************
//propert map handling
//********************

/**
 * @brief Get vertex property information
 *
 * This traits struct is used to get property information regarding AccessGraph vertices. It
 * allows access to the local descriptor, the mpl property sequence and the position, at which the
 * fusion property sequence is stored inside the bundle. It's used to allow generic property
 * selection in combination with @ref property_selector
 **/
template<typename Graph>
struct vertex_selector {
    typedef typename dcm::graph::LocalVertex key_type;
    typedef typename Graph::vertex_properties sequence_type;
    typedef mpl::int_<1> property_distance;
};

/** @brief Get edge property information
 *
 * This traits struct is used to get property information regarding AccessGraph edges. It
 * allows access to the local descriptor, the mpl property sequence and the position, at which the
 * fusion property sequence is stored inside the bundle. It's used to allow generic property
 * selection in combination with @ref property_selector
 **/
template<typename Graph>
struct edge_selector {
    typedef dcm::graph::LocalEdge key_type;
    typedef typename Graph::edge_properties sequence_type;
    typedef mpl::int_<0> property_distance;
};

/**
 * @brief Select property information trait depending on property type
 *
 * Allows generic access to property information like descriptor or property sequence by exposing
 * a specific selector type ( @ref vertex_selector or @ref edge_selector ) dependend on the supplied
 * property.
 **/
template< typename Prop, typename Graph>
struct property_selector :  mpl::if_ <
        boost::is_same <
        typename mpl::find<typename Graph::edge_properties, Prop>::type,
        typename mpl::end<typename Graph::edge_properties>::type > ,
        vertex_selector<Graph>,
        edge_selector<Graph> >::type {};

        
/** @addtogroup Metafunctions
 * @{*/
/**
 * @brief Expose if this is a edge property
 **/
template<typename T, typename Graph>
struct is_edge_property : mpl::not_< boost::is_same <
        typename mpl::find<typename Graph::edge_properties, T>::type,
        typename mpl::end<typename Graph::edge_properties>::type > > {};
/**
 * @brief Expose if this is a global edge property
 **/
template<typename T, typename Graph>
struct is_globaledge_property : mpl::not_< boost::is_same <
        typename mpl::find<typename Graph::globaledge_properties, T>::type,
        typename mpl::end<typename Graph::globaledge_properties>::type > > {};
/**
 * @brief Expose if this is a vertex property
 **/
template<typename T, typename Graph>
struct is_vertex_property : mpl::not_< boost::is_same <
        typename mpl::find<typename Graph::vertvertexex_properties, T>::type,
        typename mpl::end<typename Graph::vertex_properties>::type > > {};
/**
 * @brief Expose if this is a cluster property
 **/
template<typename T, typename Graph>
struct is_cluster_property : mpl::not_< boost::is_same <
        typename mpl::find<typename Graph::cluster_properties, T>::type,
        typename mpl::end<typename Graph::cluster_properties>::type > > {};

/**@}*/


/**
 * @brief Adapter to use dcm vertex and edge properties as boost property_maps in bgl algorithms
 *
 * Boost graph algorithms use property maps as generic storage for process data. In most cases the user
 * needs to provide these maps. The simplest way is to create them on the stack right before their
 * usage. If, however, the stored information is of use and one wants to store it permanently, this way
 * is not practical. Therefor vertex and edge properties were introduced, they allow to store arbitrary
 * information at their entity. To use this in combination with boost graph algorithms, this class can
 * be used to expose vertex and edge properties as propertie maps to the boost algorithms. All process
 * information is then stored permanently at the relevant position.
 **/
template <typename Property, typename Graph>
class property_map  {

public:
    //expose boost property map interface types
    typedef typename dcm::graph::property_selector<Property, Graph>::key_type key_type;
    typedef typename Property::type value_type;
    typedef typename Property::type&  reference;
    typedef boost::lvalue_property_map_tag category;

    //expose cutom types for easy access
    typedef Property property;
    typedef typename dcm::graph::property_selector<Property, Graph>::sequence_type sequence;
    typedef typename dcm::graph::property_selector<Property, Graph>::property_distance distance;

    /**
     * @brief Links property map with the AccessGraph which shall be accessed
     *
     * As boost graph algorithms work with local descriptors, the property map needs to know in which
     * graph they are valid. this graph has to be passed to the map. Of course this has to be the one
     * on which the algorithm is used on
     *
     * @param g shared ptr of the cluster graph on which the algorithm is used
     **/
    property_map(std::shared_ptr<Graph> g)
        : m_graph(g) { }

    std::shared_ptr<Graph> m_graph;
};


//***************************************
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
struct get_accessible_helper {

    typedef typename prop::type base_type;
    typedef base_type& result_type;
 
    template<typename bundle>
    result_type operator()(bundle& p) {
        p.template markPropertyChanged<prop>();
        return p.template getPropertyAccessible<prop>();
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

//***********************
//Function implementation
//***********************

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop, template<class, class, class, class, class> class graph_base>
template<typename T>
typename AccessGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop, graph_base>::global_extractor::result_type
AccessGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop, graph_base>::global_extractor::operator()(T& bundle) const {
    return bundle.getProperty<EdgeProperty>();
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop, template<class, class, class, class, class> class graph_base>
AccessGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop, graph_base>::global_vertex_extractor::global_vertex_extractor(AccessGraph& g) : graph(g) {};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop, template<class, class, class, class, class> class graph_base>
typename AccessGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop, graph_base>::global_vertex_extractor::result_type
AccessGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop, graph_base>::global_vertex_extractor::operator()(LocalVertex& v) const {
    return graph.getGlobalVertex(v);
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop, template<class, class, class, class, class> class graph_base>
template<typename T>
bool AccessGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop, graph_base>::operator== (const T& other) const {
    return this == &other;
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop, template<class, class, class, class, class> class graph_base>
template<typename T>
bool AccessGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop, graph_base>::operator!= (const T& other) const {
    return !(this == &other);
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop, template<class, class, class, class, class> class graph_base>
std::pair<typename AccessGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop, graph_base>::global_vertex_iterator, typename AccessGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop, graph_base>::global_vertex_iterator>
AccessGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop, graph_base>::globalVertices() {
    std::pair<local_vertex_iterator, local_vertex_iterator> res = boost::vertices(m_graph);
    global_vertex_iterator begin = boost::make_transform_iterator(res.first, global_vertex_extractor(m_graph));
    global_vertex_iterator end   = boost::make_transform_iterator(res.second, global_vertex_extractor(m_graph));

    return std::pair<global_vertex_iterator, global_vertex_iterator> (begin, end);
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop, template<class, class, class, class, class> class graph_base>
std::pair<LocalEdge, bool>
AccessGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop, graph_base>::edge(LocalVertex source, LocalVertex target) {
    return boost::edge(source, target, m_graph);
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop, template<class, class, class, class, class> class graph_base>
std::pair<LocalEdge, bool>
AccessGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop, graph_base>::edge(GlobalVertex source, GlobalVertex target) {
    
    std::pair<LocalVertex, bool> r1 = getLocalVertex(source);
    std::pair<LocalVertex, bool> r2 = getLocalVertex(target);
      
    if(!r1.second || !r2.second)
        return std::make_pair(LocalEdge(), false);
        
    return edge(r1.first, r2.first);
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop, template<class, class, class, class, class> class graph_base>
std::pair<typename AccessGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop, graph_base>::global_edge_iterator, typename AccessGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop, graph_base>::global_edge_iterator>
AccessGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop, graph_base>::getGlobalEdges(LocalEdge e) {

    //note that we can the global edges via the accessible property function as we are sure we do not change it.
    //changes to the global edges itself are tracked by them.
    auto& vec = m_graph[e].template getPropertyAccessible<GEdgeProperty>();
    
    global_edge_iterator begin = boost::make_transform_iterator(vec.begin(), global_extractor());
    global_edge_iterator end   = boost::make_transform_iterator(vec.end(), global_extractor());

    return std::pair<global_edge_iterator, global_edge_iterator> (begin, end);
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop, template<class, class, class, class, class> class graph_base>
int AccessGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop, graph_base>::getGlobalEdgeCount(LocalEdge e) {

    return m_graph[e].template getProperty<GEdgeProperty>().size();
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop, template<class, class, class, class, class> class graph_base>
std::pair<LocalEdge, bool>
AccessGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop, graph_base>::getLocalEdge(GlobalEdge e) {
    
    auto iter = boost::edges(m_graph);
    for(; iter.first != iter.second;++iter.first) {
    
        auto giter = getGlobalEdges(*iter.first);
        for(; giter.first != giter.second; ++giter.first) {
            if(*giter.first == e) 
                return std::make_pair(*iter.first, true);
        }
    }
    return std::make_pair(LocalEdge(), false);
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop, template<class, class, class, class, class> class graph_base>
GlobalVertex
AccessGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop, graph_base>::getGlobalVertex(LocalVertex v) const {
    return m_graph[v].template getProperty<VertexProperty>();
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop, template<class, class, class, class, class> class graph_base>
std::pair<LocalVertex, bool>
AccessGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop, graph_base>::getLocalVertex(GlobalVertex vertex) {
    
    auto  it = boost::vertices(m_graph);

    for(; it.first != it.second; it.first++) {
        if(vertex == m_graph[*it.first].template getProperty<VertexProperty>())
            return std::make_pair(*it.first, true);
    }
    return std::make_pair(LocalVertex(), false);
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop, template<class, class, class, class, class> class graph_base>
template<typename property, typename key>
const typename property::type&
AccessGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop, graph_base>::getProperty(key k) {
    return apply_to_bundle(k, get_property_helper<property, AccessGraph>());
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop, template<class, class, class, class, class> class graph_base>
template<typename property, typename key>
typename property::type&
AccessGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop, graph_base>::getPropertyAccessible(key k) {
    return apply_to_bundle(k, get_accessible_helper<property, AccessGraph>());
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop, template<class, class, class, class, class> class graph_base>
template<typename property, typename key>
void AccessGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop, graph_base>::setProperty(key k, const typename property::type& val) {
    apply_to_bundle(k, set_property_helper<property, AccessGraph>(val));
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop, template<class, class, class, class, class> class graph_base>
template<typename property, typename key>
bool AccessGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop, graph_base>::isPropertyChanged(key k) {
    return apply_to_bundle(k, change_property_helper<property, AccessGraph>());
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop, template<class, class, class, class, class> class graph_base>
template<typename property, typename key>
void AccessGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop, graph_base>::acknowledgePropertyChange(key k) {
    apply_to_bundle(k, acknowledge_property_helper<property, AccessGraph>());
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop, template<class, class, class, class, class> class graph_base>
template<typename key>
bool AccessGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop, graph_base>::hasPropertyChanges(key k) {
    return apply_to_bundle(k, changes_property_helper<AccessGraph, key>());
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop, template<class, class, class, class, class> class graph_base>
template<typename key>
void AccessGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop, graph_base>::acknowledgePropertyChanges(key k) {
    apply_to_bundle(k, ack_changes_property_helper<AccessGraph, key>());
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop, template<class, class, class, class, class> class graph_base>
bool AccessGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop, graph_base>::edgeChanged(LocalEdge e) {
  
    if(m_graph[e].hasPropertyChanges())
        return true;
    
    auto& vec = m_graph[e].template getProperty<GEdgeProperty>();
    for(const edge_bundle_single& global : vec) {
        if(global.template hasPropertyChanges())
            return true;
    }
    return false;
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop, template<class, class, class, class, class> class graph_base>
void AccessGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop, graph_base>::acknowledgeEdgeChanges(LocalEdge e) {
  
    m_graph[e].acknowledgePropertyChanges();
        
    auto& vec = m_graph[e].template getPropertyAccessible<GEdgeProperty>();
    for(edge_bundle_single& global : vec)
        global.acknowledgePropertyChanges();
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop, template<class, class, class, class, class> class graph_base>
void AccessGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop, graph_base>::initIndexMaps() {

    //just iterate over all edges and vertices and give them all a unique index
    std::pair<local_vertex_iterator, local_vertex_iterator>  vit = boost::vertices(m_graph);

    for(int c = 0; vit.first != vit.second; ++vit.first, ++c)
        setProperty<Index>(*vit.first, c);

    std::pair<local_edge_iterator, local_edge_iterator>  eit = boost::edges(m_graph);

    for(int c = 0; eit.first != eit.second; ++eit.first, ++c)
        setProperty<Index>(*eit.first, c);
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop, template<class, class, class, class, class> class graph_base>
template<typename functor>
typename functor::result_type
AccessGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop, graph_base>::apply_to_bundle(LocalVertex k, functor f) {
    return f(m_graph[k]);
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop, template<class, class, class, class, class> class graph_base>
template<typename functor>
typename functor::result_type
AccessGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop, graph_base>::apply_to_bundle(LocalEdge k, functor f) {
    return f(m_graph[k]);
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop, template<class, class, class, class, class> class graph_base>
template<typename functor>
typename functor::result_type
AccessGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop, graph_base>::apply_to_bundle(GlobalVertex k, functor f) {

    //check all vertices if they are the id
    std::pair<local_vertex_iterator, local_vertex_iterator>  it = boost::vertices(m_graph);

    for(; it.first != it.second; it.first++) {
        vertex_bundle& p = m_graph[*it.first];

        if(k == p.template getProperty<VertexProperty>())
            return f(p);
    }

    dcm_assert(false);
        //TODO: Throw (propeties return reference, but cant init a reference temporarily)
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop, template<class, class, class, class, class> class graph_base>
template<typename functor>
typename functor::result_type
AccessGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop, graph_base>::apply_to_bundle(GlobalEdge k, functor f) {

    auto iter = boost::edges(m_graph);
    for(; iter.first != iter.second;++iter.first) {
    
        //search the global one in the local edge. Note that we can easily use the accessible version 
        //as we dont change the vector. If the functor changes the individual properties we this will be tracked
        //by them individual
        auto& vec = m_graph[*iter.first].template getPropertyAccessible<GEdgeProperty>();
        edge_single_iterator it;

        for(it = vec.begin(); it != vec.end(); it++) {
            if((*it).template getProperty<EdgeProperty>() == k)
                return f(*it);
        };
    }

    dcm_assert(false)
    //TODO: throw, as there has to be a global edge and the following return statement is illegal
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop, template<class, class, class, class, class> class graph_base>
template<typename Predicate, typename Iterator>
std::pair<boost::filter_iterator<Predicate, Iterator>, boost::filter_iterator<Predicate, Iterator> >
AccessGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop, graph_base>::filterRange(std::pair<Iterator, Iterator> range) {

    Predicate p = Predicate(sp_base::shared_from_this());
    return std::make_pair(boost::make_filter_iterator(p, range.first, range.second),
                          boost::make_filter_iterator(p, range.second, range.second));
};
    
} //namespace graph
} //namespace dcm


namespace boost {
//access the propertymap needs to be boost visable
template<typename P, typename G>
typename dcm::graph::property_map<P, G>::value_type    get(const dcm::graph::property_map<P, G>& map,
        typename dcm::graph::property_map<P, G>::key_type key)
{

    typedef dcm::graph::property_map<P, G> map_t;
    return  map.m_graph->template getProperty<P>(key);
};

template <typename P, typename G>
void  put(const dcm::graph::property_map<P, G>& map,
          typename dcm::graph::property_map<P, G>::key_type key,
          const typename dcm::graph::property_map<P, G>::value_type& value)
{

    typedef dcm::graph::property_map<P, G> map_t;
    map.m_graph->template setProperty<P>(key, value);
};


template <typename P, typename G>
typename dcm::graph::property_map<P, G>::reference at(const dcm::graph::property_map<P, G>& map,
        typename dcm::graph::property_map<P, G>::key_type key)
{
    typedef dcm::graph::property_map<P, G> map_t;
    return map.m_graph->template getPropertyAccessible<P>(key);
}
}


#endif // ACCESSGRAPH_HPP





