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

#include "signal.hpp"
#include "property.hpp"


namespace dcm {

typedef int ObjectTypeID;    
//few standart signal names
struct remove {};

namespace details {

#define EMIT_OBJECT_SIGNAL_CALL_DEC(z, n, data) \
    template<typename SigMap> \
    template < \
    typename S  \
    BOOST_PP_ENUM_TRAILING_PARAMS(n, typename Arg) \
    > \
    typename boost::enable_if<mpl::has_key<SigMap, S>, void>::Type \
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

    /** @defgroup Objects Objects
 *
 * @brief Concept and functionality of the dcm objects
 *
 *
 **/
template<typename Final>
struct Object {
  
    Object();
    
    //functions for accessing stacked properties. we basicly mirror the PropertyOwner interface
    //so that it looks to the user like we actually are a propertyowner but then access the propertyowner 
    //which holds the property for query
    template<typename Prop>
    const typename Prop::type& getProperty();    
    template<typename Prop>
    typename boost::disable_if<details::has_change_tracking<Prop> >::type setProperty(const typename Prop::type& value);
    template<typename Prop>
    typename boost::enable_if<details::has_change_tracking<Prop> >::type setProperty(const typename Prop::type& value);
    
    template<typename S>
    Connection connectSignal(typename mpl::at<SigMap, S>::type function);
    template<typename S>
    void disconnectSignal(Connection c);
    
    //object identification 
    const ObjectTypeID getTypeID();
    template<typename Object>
    bool isObject();
    
protected:
    //with no vararg templates before c++11 we need preprocessor to create the overloads of emit signal we need
    BOOST_PP_REPEAT(5, EMIT_SIGNAL_CALL_DEF, ~)
    
    const ObjectTypeID m_id;
};

template<typename Final>
template<typename Prop>
const typename Prop::type& Object<Final>::getProperty() {
    
       typedef mpl::at_key<Prop, typename Final::PropertyMap>::type Object;
       
       if( Final::ObjectTypeID<Object>::value == m_id ) 
           return static_cast<Object*>(this)->m_properties->getProperty<Prop>();
       //else 
           //TODO: throw
};

}; //details
}; //dcm

#endif //DCM_SIGNAL_H


