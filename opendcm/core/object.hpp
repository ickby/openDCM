/*
    openGCM, geometric constraint manager
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

#ifndef GCM_OBJECT_H
#define GCM_OBJECT_H

#include <iostream>
#include <list>

#include <boost/mpl/at.hpp>
#include <boost/mpl/if.hpp>
#include <boost/mpl/void.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/mpl/transform.hpp>

#include <boost/fusion/include/as_vector.hpp>
#include <boost/fusion/include/mpl.hpp>
#include <boost/fusion/include/at.hpp>

#include <boost/preprocessor.hpp>
#include <boost/preprocessor/repetition/repeat.hpp>
#include <boost/preprocessor/cat.hpp>
#include <boost/preprocessor/repetition/enum.hpp>
#include <boost/preprocessor/repetition/enum_params.hpp>
#include <boost/preprocessor/repetition/enum_trailing_params.hpp>
#include <boost/preprocessor/repetition/enum_binary_params.hpp>

#include <boost/enable_shared_from_this.hpp>
#include <boost/any.hpp>

#include "property.hpp"


namespace mpl = boost::mpl;
namespace fusion = boost::fusion;

namespace dcm {

namespace details {
template<typename T>
struct map_key {
    typedef typename T::first type;
};
template<typename T>
struct map_val {
    typedef typename T::second type;
};

}

typedef boost::any Connection;

/**
 * @brief Base class for all object types
 *
 * This class add's property and signal capabilitys to all deriving classes. For properties it is tigthly
 * integrated with the system class: It searches systems property list for the derived class as specified by
 * the second template parameter and makes it accessible via appopriate functions. Signals are speciefied by a
 * mpl::map with signal name type as key and a boost::function as values.
 *
 * \tparam Sys class of type System of which the properties are taken
 * \tparam Obj the type of the derived object
 * \tparam Sig a mpl::map specifing the object's signals by (type -  boost::function) pairs
 **/
template<typename Sys, typename Obj, typename Sig>
struct Object : public boost::enable_shared_from_this<Obj> {

    Object(Sys& system) : m_system(system) {};

    /**
      * @brief Access properties
      *
      * Returns a reference to the propertys actual value. The property type has to be registerd to the
      * System type which was given as template parameter to this object.
      * @tparam Prop property type which should be accessed
      * @return Prop::type& a reference to the properties actual value.
      **/
    template<typename Prop>
    typename Prop::type& getProperty() {
        typedef typename mpl::find<Sequence, Prop>::type iterator;
        typedef typename mpl::distance<typename mpl::begin<Sequence>::type, iterator>::type distance;
        return fusion::at<distance>(m_properties);
    };

    /**
       * @brief Set properties
       *
       * Set'S the value of a specified property. The property type has to be registerd to the
       * System type which was given as template parameter to this object. Note that setProperty(value)
       * is equivalent to getProperty() = value.
       * @tparam Prop property type which should be setProperty
       * @param value value of type Prop::type which should be set in this object
       **/
    template<typename Prop>
    void setProperty(typename Prop::type value) {
        typedef typename mpl::find<Sequence, Prop>::type iterator;
        typedef typename mpl::distance<typename mpl::begin<Sequence>::type, iterator>::type distance;
        fusion::at<distance>(m_properties) = value;
    };

    /**
     * @brief Connects a slot to a specified signal.
     *
     * Slots are boost::functions which get called when the signal is emitted. Any valid boost::function
     * which ressembles the signal tyes signature can be registert. It is important that the signal type
     * was registerd to this object on creation by the appropriate template parameter.
     *
     * @tparam S the signal which should be intercepted
     * @param function boost::function which resembles the signal type's signature
     * @return void
     **/
    template<typename S>
    Connection connectSignal(typename mpl::at<Sig, S>::type function) {
        typedef typename mpl::find<sig_name, S>::type iterator;
        typedef typename mpl::distance<typename mpl::begin<sig_name>::type, iterator>::type distance;
        typedef typename fusion::result_of::value_at<Signals, distance>::type list_type;
        list_type& list = fusion::at<distance>(m_signals);
        return list.insert(list.begin(),function);
    };

    /**
    * @brief Disconnects a slot for a specific signal.
    *
    * Disconnects a slot so that it dosn't get called at signal emittion. It's important to
    * disconnect the slot by the same boost:function it was connected with.
    *
    * @tparam S the signal type of interest
    * @param function boost::function with which the slot was connected
    * @return void
    **/
    template<typename S>
    void disconnectSignal(Connection c) {
        typedef typename mpl::find<sig_name, S>::type iterator;
        typedef typename mpl::distance<typename mpl::begin<sig_name>::type, iterator>::type distance;

        typedef typename fusion::result_of::value_at<Signals, distance>::type list_type;
        list_type& list = fusion::at<distance>(m_signals);
        list.erase(boost::any_cast<typename list_type::iterator>(c));
    };

protected:

    /*properties
     * search the property map of the system class and get the mpl::vector of properties for the
     * derived type. It's imortant to not store the properties but their types. These types are
     * stored and accessed as fusion vector.
     * */
    typedef typename mpl::at<typename Sys::object_properties, Obj>::type Mapped;
    typedef typename mpl::if_< boost::is_same<Mapped, mpl::void_ >, mpl::vector<>, Mapped>::type Sequence;
    typedef typename mpl::transform<Sequence, details::property_type<mpl::_1> >::type Typesequence;
    typedef typename fusion::result_of::as_vector<Typesequence>::type Properties;

    /*signal handling
     * extract all signal types to allow index search (inex search on signal functions would fail as same
     * signatures are supported for multiple signals). Create std::vectors to allow multiple slots per signal
     * and store these vectors in a fusion::vector for easy access.
     * */
    typedef typename mpl::fold< Sig, mpl::vector<>,
            mpl::push_back<mpl::_1, details::map_key<mpl::_2> > >::type sig_name;
    typedef typename mpl::fold< Sig, mpl::vector<>,
            mpl::push_back<mpl::_1, details::map_val<mpl::_2> > >::type sig_functions;
    typedef typename mpl::fold< sig_functions, mpl::vector<>,
            mpl::push_back<mpl::_1, std::list<mpl::_2> > >::type sig_vectors;
    typedef typename fusion::result_of::as_vector<sig_vectors>::type Signals;

    Sys& m_system;
    Properties m_properties;
    Signals m_signals;

    //with no vararg templates before c++11 we need preprocessor to create the overloads of emit signal we need
#define EMIT_ARGUMENTS(z, n, data) \
    BOOST_PP_CAT(data, n)

#define EMIT_CALL(z, n, data) \
    template < \
    typename S  \
    BOOST_PP_ENUM_TRAILING_PARAMS(n, typename Arg) \
    > \
    void emitSignal( \
                     BOOST_PP_ENUM_BINARY_PARAMS(n, Arg, const& arg) \
                   ) \
    { \
        typedef typename mpl::find<sig_name, S>::type iterator; \
        typedef typename mpl::distance<typename mpl::begin<sig_name>::type, iterator>::type distance; \
        typedef typename fusion::result_of::value_at<Signals, distance>::type list_type; \
        list_type& list = fusion::at<distance>(m_signals); \
        for (typename list_type::iterator it=list.begin(); it != list.end(); it++) \
            (*it)(BOOST_PP_ENUM(n, EMIT_ARGUMENTS, arg)); \
    };

    BOOST_PP_REPEAT(5, EMIT_CALL, ~)



};

}

#endif //GCM_OBJECT_H

