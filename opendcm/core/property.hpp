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

#ifndef GCM_PROPERTY_H
#define GCM_PROPERTY_H

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/buffer_concepts.hpp>
#include <boost/graph/properties.hpp>
#include <boost/fusion/sequence.hpp>
#include <boost/fusion/container/vector.hpp>

#include <boost/mpl/find.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/property_map/property_map.hpp>

namespace mpl = boost::mpl;
namespace fusion = boost::fusion;

namespace dcm {
  
/** @addtogroup Core
 * @{ */

/** @defgroup Property Properties
 * 
 * @brief Concept and handling of properties for generic data storage and type extension.
 * 
 * Properties are a basic building block of the dcm and fullfill two essential tasks: First, build a 
 * infrastructure for storing data in any kind of object and second, make this data universally accessible.
 * Universally accessible means in this context, that it shall be possible to retrieve data only by 
 * knowing some kind of identifier. No data type specific get/set functions or access to class members
 * should be needed to access the stored values. The usage of identifiers allows to design interfaces
 * for properties in a type interchangable way. Therefore no restrictions are imposed on the interface,
 * no matter what or how much data is stored.
 * 
 * The connection of data type and identifier is achieved through the property structs, which all follow
 * the same concept: Identifier is the struct type, the stored data is exposed as 'type' typedef. The data 
 * type can be every c++ type (including classes and structs) which is default constructable. They don't need
 * to be assignable or copyable by default, thats only nesseccary if you want to change the hole stored 
 * object by assigning or set-methods. If not, the data object can be uncopyable and it should be used by 
 * retrieving it's reference with get-methods.
 * 
 * Propertys are further designed to fit in the concept of compile-time modularisation. To allow the extension 
 * of all data-holding entitys with new data types, propertys store their own purpose. Thats 
 * done by extending the property struct with a second typedef which is named kind and which specifies of which
 * kind the property is. That means, that this typedef defines when the property shall be used and for which 
 * context it is designed for. Dependend on the propertys kind, it will be added to diffrent places inside the dcm.
 * A property of kind @ref vertex_property will added to vertices, a property of kind @ref object_property to all
 * objects and so on. A property implementation for storing integers at a graph edge with the identifier
 * 'test'property' may look like that:
 * @code 
 * struct test_property {
 * 	typedef int type;
 * 	typedef edge_property kind;
 * }
 * @endcode
 * 
 * @{ */

/**
 * @brief Identifier for vertex properties
 * 
 * This is a identifier structure for vertex properties. Every property with this struct as 'kind' type
 * will be added to all vertices of a cluster. It is accessible through global and local vertex
 * descriptors. These properties are intended for use in boost algorithms in combination with 
 * \ref property_map .
 */
struct vertex_property {};
/**
 * @brief Identifier for edge properties
 * 
 * This is a identifier structure for edge properties. Every property with this struct as 'kind' type
 * will be added to all local edges of a cluster. It is accessible through local edge
 * descriptors, or global one by getting it's holding local edge first. Note that global edges don't 
 * have properties, as the properties are intended for use inside boost graph algorithms and therefore 
 * only needed in local edges. @see property_map
 */
struct edge_property {};

/**
 * @brief Identifier for cluster properties
 * 
 * A ClusterGraph has it's own properties, and ever property with this identifier as 'kind' type will be
 * added to it. This is intended for internal dcm usage, its possible to give the abstract cluster a meaning
 * by adding special properties to it. It can be accessed by special ClusterGraph functions designed for this
 * purpose.
 **/
struct cluster_property {};

/**
 *@brief Identifier for general object properties
 *
 * Aproperty with this struct as 'kind' type will be added to all existing objects, no matter of individual
 * type. Use this only for general, sharable properties. To add a property to a single object, use it's 
 * type as 'kind'.
 **/
struct object_property {};

namespace details {

/** @addtogroup Metafunctions
 * @{
 * @brief Get vertex property information 
 * 
 * This traits struct is used to get property information regarding ClusterGraph vertices. It
 * allows access to the local descriptor, the mpl property sequence and the position, at which the 
 * fusion property sequence is stored inside the bundle. It's used to allow generic property 
 * selection in combination with @ref property_selector
 **/
template<typename Graph>
struct vertex_selector {
    typedef typename boost::graph_traits<Graph>::vertex_descriptor key_type;
    typedef typename Graph::vertex_properties sequence_type;
    typedef mpl::int_<1> property_distance;
};

/** @brief Get edge property information 
 * 
 * This traits struct is used to get property information regarding ClusterGraph edges. It
 * allows access to the local descriptor, the mpl property sequence and the position, at which the 
 * fusion property sequence is stored inside the bundle. It's used to allow generic property 
 * selection in combination with @ref property_selector
 **/
template<typename Graph>
struct edge_selector {
    typedef typename boost::graph_traits<Graph>::edge_descriptor key_type;
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
template< typename Kind, typename Graph>
struct property_selector : public mpl::if_<boost::is_same<Kind, vertex_property>,
        vertex_selector<Graph>, edge_selector<Graph> >::type {};

/**
 * @brief Metafunction to expose the property storage type
 **/
template<typename T>
struct property_type {
    typedef typename T::type type;
};
/**
 * @brief Metafunction to expose which kid of property this is
 **/
template<typename T>
struct property_kind {
    typedef typename T::kind type;
};

/**
 * @brief Property vector to a fusion sequence of the propety storage types
 * 
 * Properties are passed around as mpl sequences, mostly vectors. To store actual values, they need to
 * be transformed into fusion sequences. However, only the storage type needs to be in the vector, not
 * the 'kind' information (as this is only used as meta information for type creation). This struct 
 * exposes a fusion vector which can hold all property storage types.
 **/
template<typename T>
struct pts { //property type sequence
    typedef typename mpl::transform<T, details::property_type<mpl::_1> >::type ptv;
    typedef typename fusion::result_of::as_vector< ptv >::type type;
};
/**@}*/
}

/** @addtogroup Metafunctions
 * @{
/**
 * @brief Expose if this is a edge property
 **/
template<typename T>
struct is_edge_property : boost::is_same<typename T::kind,edge_property> {};
/**
 * @brief Expose if this is a vertex property
 **/
template<typename T>
struct is_vertex_property : boost::is_same<typename T::kind,vertex_property> {};
/**
 * @brief Expose if this is a cluster property
 **/
template<typename T>
struct is_cluster_property : boost::is_same<typename T::kind,cluster_property> {};
/**
 * @brief Expose if this is a general object property
 **/
template<typename T>
struct is_object_property : boost::is_same<typename T::kind,object_property> {};
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
    typedef typename dcm::details::property_selector<typename Property::kind, Graph>::key_type key_type;
    typedef typename Property::type value_type;
    typedef typename Property::type&  reference;
    typedef boost::lvalue_property_map_tag category;

    //expose cutom types for easy access
    typedef Property property;
    typedef typename dcm::details::property_selector<typename Property::kind, Graph>::sequence_type sequence;
    typedef typename dcm::details::property_selector<typename Property::kind, Graph>::property_distance distance;

    /**
     * @brief Links property map with the ClusterGraph which shall be accessed
     * 
     * As boost graph algorithms work with local descriptors, the property map needs to know in which
     * graph they are valid. this graph has to be passed to the map. Of course this has to be the one
     * on which the algorithm is used on
     *
     * @param g shared ptr of the cluster graph on which the algorithm is used 
     **/
    property_map(boost::shared_ptr<Graph> g)
        : m_graph(g) { }

    boost::shared_ptr<Graph> m_graph;
};


//now create some standart properties
//***********************************

/**
 * @brief Dummy property
 **/
struct empty_prop {
  typedef int kind;
  typedef int type;
};
//type of a graph cluster
/**
 * @brief Add a type to clusters
 * 
 * Allows to specify special types to ClusterGraphs and make a it possibe to distuingish between
 * diffrent purposes. The cluster types need to be int.
 **/
struct type_prop {
    //states the type of a cluster
    typedef cluster_property kind;
    typedef int type;
};
//cluster in graph changed?
/**
 * @brief Was the cluster changed?
 * 
 * Adds a boolean to the cluster which indicates if the cluster was changedsince the last
 * processing. It should be set to true if vertices and edges were added or removed, Subclusters
 * created or deleted and so on.
 **/
struct changed_prop {
    typedef cluster_property kind;
    typedef bool type;
};
/**
 * @brief Add an index to vertices
 * 
 * Most boost graph algorithms need a index for vertices, ranging from 0 to vertex count. As this can
 * be useful for many other things it is added as vertex property.
 **/
struct vertex_index_prop {
    typedef vertex_property kind;
    typedef int type;
};
/**
 * @brief Add an index to edges
 * 
 * Most boost graph algorithms need a index for edges, ranging from 0 to edge count. As this can
 * be useful for many other things it is added as edge property.
 **/
struct edge_index_prop {
    typedef edge_property kind;
    typedef int type;
};
/**
 * @brief Add an ID to objects
 * 
 * It may be wanted to add identification markers to objects, this property can be used for that. It 
 * is special, as it takes its storage type as template parameter. This property can therefore not be
 * directly accessed with this struct, the template parameter has to be known.
 * @tparam T the identifier type
 **/
template<typename T>
struct id_prop {
    typedef object_property kind;
    typedef T type;
};

/**@}*/ 
/**{*/
}

namespace boost {
//access the propertymap needs to be boost visable
template<typename P, typename G>
typename dcm::property_map<P,G>::value_type	get(const dcm::property_map<P,G>& map,
        typename dcm::property_map<P,G>::key_type key)  {

    typedef dcm::property_map<P,G> map_t;
    typedef typename mpl::find<typename map_t::sequence, typename map_t::property>::type iterator;
    typedef typename mpl::distance<typename mpl::begin<typename map_t::sequence>::type, iterator>::type distance;

    return  fusion::at<distance>(fusion::at<typename dcm::property_map<P,G>::distance>(map.m_graph->operator[](key)));
};

template <typename P, typename G>
void  put(const dcm::property_map<P,G>& map,
          typename dcm::property_map<P,G>::key_type key,
          const typename dcm::property_map<P,G>::value_type& value)  {

    typedef dcm::property_map<P,G> map_t;
    typedef typename mpl::find<typename map_t::sequence, typename map_t::property>::type iterator;
    typedef typename mpl::distance<typename mpl::begin<typename map_t::sequence>::type, iterator>::type distance;
    fusion::at<distance>(fusion::at<typename dcm::property_map<P,G>::distance>(map.m_graph->operator[](key))) = value;
};


template <typename P, typename G>
typename dcm::property_map<P,G>::reference at(const dcm::property_map<P,G>& map,
        typename dcm::property_map<P,G>::key_type key) {
    typedef dcm::property_map<P,G> map_t;
    typedef typename mpl::find<typename map_t::sequence, typename map_t::property>::type iterator;
    typedef typename mpl::distance<typename mpl::begin<typename map_t::sequence>::type, iterator>::type distance;
    return fusion::at<distance>(fusion::at<typename dcm::property_map<P,G>::distance>(map.m_graph->operator[](key)));
}
}

#endif //GCM_PROPERTY_H
