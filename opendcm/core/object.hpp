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
#include <boost/mpl/range_c.hpp>

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
struct Builder;
}

namespace details {  

//if a object was derived from a PropertyOwner than the Properties vector is not empty
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
    virtual void preprocessVertex(std::shared_ptr<graph::AccessGraphBase>, graph::LocalVertex, graph::GlobalVertex) {};
    virtual void preprocessEdge(std::shared_ptr<graph::AccessGraphBase>, graph::GlobalEdge) {};
    
    virtual void postprocessVertex(std::shared_ptr<graph::AccessGraphBase>, graph::LocalVertex, graph::GlobalVertex) {};
    virtual void postprocessEdge(std::shared_ptr<graph::AccessGraphBase>, graph::GlobalEdge) {};
    
    friend struct solver::Builder;
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
/*
    typedef mpl::vector0<> Properties;
    //functions for accessing properties without the need to transform the object manual
    template<typename Prop>
    const typename Prop::type& getProperty();
    template<typename Prop>
    void setProperty(const typename Prop::type& value);
 */
    /*template<typename S>
    Connection connectSignal(typename mpl::at<SigMap, S>::type function);
    template<typename S>
    void disconnectSignal(Connection c);*/

    //object identification
    const ObjectTypeID getTypeID();
    template<typename Object>
    bool isType();

protected:    
    const ObjectTypeID m_id;
};

template<typename Final>
Object<Final>::Object(int ID) : m_id(ID) {

};
template<typename T>
void pretty(T t) {
    std::cout<<__PRETTY_FUNCTION__<<std::endl;
};

template<typename List, typename Prop>
struct objGetProp {
  
    objGetProp(void* obj, int id) : m_object(obj) {};
    
    template<typename T>
    void operator()(const T& t) const {
        if(T::value == id)
            m_store = &static_cast<typename mpl::at<List, T>::type>(m_object)->template getProperty<Prop>();
    }
    
    const typename Prop::type& getResult() {return m_store;};
    
private:
    typename Prop::type* m_store = nullptr;
    void* m_object = nullptr;
    int id;
};
/*
template<typename Final>
template<typename Prop>
const typename Prop::type& Object<Final>::getProperty() {
   
    //iterate over all base objects to find the one we are
    typedef mpl::range_c<int,0,
            mpl::size<typename Final::ObjectList>::value> Range;
    auto func = objGetProp<typename Final::ObjectList, Prop>(this, getTypeID());
    mpl::for_each<Range>(func);
    return func->getResult();
};

template<typename Final>
template<typename Prop>
void Object<Final>::setProperty(const typename Prop::type& value) {
    
    typedef typename Final::template objectByProperty<Prop>::type Object;
    //filteres list of all types inheriting from Object
    typedef typename mpl::remove_if<
    typename Final::ObjectList,
             mpl::not_<boost::is_base_of<Object, mpl::_> > >::type Filtered;
            
    if(isOneOfTypes<typename Final::ObjectList, Filtered>())
        static_cast<Object*>(this)->m_properties.template setProperty<Prop>(value);
    else
        throw property_error() <<  boost::errinfo_errno(3) << error_message("property does not exist in this object");
};
*/
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


