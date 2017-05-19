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

#ifndef DCM_OBJECT_H
#define DCM_OBJECT_H

#include <memory>

#include "signal.hpp"
#include "property.hpp"

#include <boost/mpl/equal_to.hpp>
#include <boost/mpl/remove_if.hpp>
#include <boost/preprocessor.hpp>
#include <boost/preprocessor/repetition/repeat.hpp>
#include <boost/preprocessor/tuple/to_seq.hpp>
#include <boost/preprocessor/tuple/size.hpp>
#include <boost/preprocessor/seq/for_each.hpp>
#include <boost/preprocessor/seq/enum.hpp>
#include <boost/preprocessor/cat.hpp>

#ifndef DCM_MAX_OBJECTS
#define DCM_MAX_OBJECTS 10
#endif

namespace mpl = boost::mpl;

namespace dcm {
    
typedef int ObjectTypeID;
//few standart signal names
struct remove {};

//we need this forward declaration to allow a friend statement later on 
namespace solver {
template<typename Kernel> struct Builder;
}

namespace details {

#define EMIT_OBJECT_SIGNAL_CALL_DEC(z, n, data) \
    template<typename SigMap> \
    template < \
    typename S  \
    BOOST_PP_ENUM_TRAILING_PARAMS(n, typename Arg) \
    > \
    typename boost::enable_if<mpl::has_key<SigMap, S>, void>::type \
    SignalOwner<SigMap>::emitSignal( \
            BOOST_PP_ENUM_BINARY_PARAMS(n, Arg, const& arg) \
                                              ) \
    { \
        typedef typename mpl::find<sig_name, S>::type iterator; \
        typedef typename mpl::distance<typename mpl::begin<sig_name>::type, iterator>::type distance; \
        typedef typename fusion::result_of::value_at<Signals, distance>::type map_type; \
        map_type& map = fusion::at<distance>(m_signals); \
        for (typename map_type::iterator it=map.begin(); it != map.end(); it++) \
            (it->second)(BOOST_PP_ENUM(n, EMIT_ARGUMENTS, arg)); \
    };

#define CHECK_TYPE(z, n, data) \
    case n :\
        return mpl::contains< Filtered, typename mpl::at<List, mpl::int_<n> >::type >::value;\
        break;
    
#define CHECK_TYPE_DEF(z, n, data) \
    template <typename List, typename Filtered> \
    typename boost::enable_if<mpl::equal_to<mpl::size<List>, mpl::int_<n> >, bool>::type \
    isOneOfTypes() \
    { \
        switch(m_id) { \
        BOOST_PP_REPEAT(n, CHECK_TYPE, ~); \
            default: \
                return false; \
        }; \
    };
    
    
//macros to ease the handling of object handling
#define PROPERTY_HANDLING_FUNCTIONS(r, data, elem) \
    public: \
        const typename elem::type& BOOST_PP_CAT(get, elem)() { \
            return m_properties.getProperty< elem >(); \
        }; \
        void BOOST_PP_CAT(set, elem)(const typename elem::type& value) { \
            m_properties.setProperty< elem >(value); \
        };
    
#define DCM_OBJECT_ADD_PROPERTIES(Base, seq) \
    BOOST_PP_SEQ_FOR_EACH(PROPERTY_HANDLING_FUNCTIONS, _, seq ) \
    typedef mpl::vector< \
        BOOST_PP_SEQ_ENUM( seq ) \
        > Properties; \
    template<typename Prop> \
    const typename std::enable_if<dcm::details::has_property<Prop, Properties>::type::value, typename Prop::type>::type& getProperty() { \
        return m_properties.getProperty<Prop>(); \
    } \
    template<typename Prop> \
    const typename std::enable_if<!dcm::details::has_property<Prop, Properties>::type::value, typename Prop::type>::type& getProperty() { \
        return Base::template getProperty<Prop>(); \
    } \
    template<typename Prop> \
    typename std::enable_if<dcm::details::has_property<Prop, Properties>::type::value>::type setProperty(const typename Prop::type& value){ \
        m_properties.setProperty<Prop>(value); \
    }; \
    template<typename Prop> \
    typename std::enable_if<!dcm::details::has_property<Prop, Properties>::type::value>::type setProperty(const typename Prop::type& value){ \
        Base::template setProperty<Prop>(value); \
    }; \
    protected: \
        dcm::details::PropertyOwner< Properties > m_properties;

template<typename Obj, typename Prop>        
struct object_has_property : public mpl::contains<typename Obj::Properties, Prop> {};
     
/**
 * @brief Obejct to be stored in the graph 
 * This object type is intended to be added to the graph for utilizing pre and postprocess operations 
 * which setup the global vertices and edges. Note that due to the graph manipulation algorithms the local 
 * descriptors can be different and also valid in different graphs for pre and post process operations.
 * It is therefore important to not store the local descriptors and graph but always use the ones passed 
 * to the functions.
 */
struct GraphObject : std::enable_shared_from_this<GraphObject> {
       
protected:
    ///preprocessing all vertices. This call happens after all subcluster of the given graph are preprocessed
    virtual void preprocessVertex(std::shared_ptr<graph::AccessGraphBase>, graph::LocalVertex, graph::GlobalVertex) {};
    ///preprocess all edges. This call happens after all subcluster of the given graph are preprocessed
    virtual void preprocessEdge(std::shared_ptr<graph::AccessGraphBase>, graph::GlobalEdge) {};
    
    ///postprocessing all vertices. This call happens after all subcluster of the given graph are preprocessed
    virtual void postprocessVertex(std::shared_ptr<graph::AccessGraphBase>, graph::LocalVertex, graph::GlobalVertex) {};
    ///postprocessing all edges. This call happens after all subcluster of the given graph are preprocessed
    virtual void postprocessEdge(std::shared_ptr<graph::AccessGraphBase>, graph::GlobalEdge) {};
    
    
    //postprocess clusters. Before postprocessing all vertices and edges within a cluster this call happens
    virtual void preprocessCluster(std::shared_ptr<graph::AccessGraphBase>,     //the cluster currently processed
                                   graph::LocalVertex,                          //the vertex within the cluster that is a subcluster
                                   std::shared_ptr<graph::AccessGraphBase>) {}; //the subcluster the vertex represents
    //postprocess clusters. Before postprocessing all vertices and edges within a cluster this call happens
    virtual void postprocessCluster(std::shared_ptr<graph::AccessGraphBase>,     //the cluster currently processed
                                    graph::LocalVertex,                          //the vertex within the cluster that is a subcluster
                                    std::shared_ptr<graph::AccessGraphBase>) {}; //the subcluster the vertex represents
    
    template<typename Kernel> friend struct solver::Builder;
};

struct GraphObjectProperty {
    typedef std::shared_ptr<GraphObject> type;
};

/** @defgroup Objects Objects
*
* @brief Concept and functionality of the dcm objects. 
* They can be stored at graphs and provide property access methods.
* 
**/
template<typename Final>
struct Object : public GraphObject {

    Object(int ID);

    //functions for accessing stacked properties. we basicly mirror the PropertyOwner interface
    //so that it looks to the user like we actually are a propertyowner but then access the propertyowner
    //which holds the property for query
    template<typename Prop>
    const typename Prop::type& getProperty();
    template<typename Prop>
    void setProperty(const typename Prop::type& value);
 
    /*template<typename S>
    Connection connectSignal(typename mpl::at<SigMap, S>::type function);
    template<typename S>
    void disconnectSignal(Connection c);*/

    //object identification
    const ObjectTypeID getTypeID();
    template<typename Object>
    bool isType();

protected:    
    //with no vararg templates before c++11 we need preprocessor to create the overloads of emit signal we need
    //BOOST_PP_REPEAT(5, EMIT_SIGNAL_CALL_DEF, ~)

    BOOST_PP_REPEAT(DCM_MAX_OBJECTS, CHECK_TYPE_DEF, ~)

    const ObjectTypeID m_id;
};

template<typename Final>
Object<Final>::Object(int ID) : m_id(ID) {

};
template<typename T>
void pretty(T t) {
    std::cout<<__PRETTY_FUNCTION__<<std::endl;
};

template<typename Final>
template<typename Prop>
const typename Prop::type& Object<Final>::getProperty() {
   
    throw property_error() <<  boost::errinfo_errno(3) << error_message("property does not exist in this object");
};

template<typename Final>
template<typename Prop>
void Object<Final>::setProperty(const typename Prop::type& value) {
    
    throw property_error() <<  boost::errinfo_errno(3) << error_message("property does not exist in this object");
};

template<typename Final>
const ObjectTypeID Object<Final>::getTypeID() {
    
    return m_id;
}

template<typename Final>
template<typename Obj>
bool Object<Final>::isType() {
    
    return Final::template objectTypeID<Obj>::ID::value == m_id;
}

}; //details
}; //dcm

#endif //DCM_SIGNAL_H


