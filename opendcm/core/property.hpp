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

#ifndef DCM_PROPERTY_H
#define DCM_PROPERTY_H

#include <boost/mpl/find.hpp>
#include <boost/mpl/void.hpp>
#include <boost/mpl/filter_view.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/mpl/for_each.hpp>
#include <boost/mpl/transform.hpp>

#include <boost/fusion/mpl.hpp>
#include <boost/fusion/include/vector.hpp>
#include <boost/fusion/include/at.hpp>
#include <boost/fusion/include/for_each.hpp>

#include <boost/phoenix/phoenix.hpp>
#include <boost/phoenix/fusion.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/exception/errinfo_errno.hpp>
#include <boost/function.hpp>

#include "defines.hpp"
#include "signal.hpp"

namespace mpl = boost::mpl;
namespace fusion = boost::fusion;

namespace dcm {

/** @addtogroup Core
 * @{ */

/** @defgroup Property Properties
 *
 * @brief Concept and handling of properties for generic data storage.
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
 * to be assignable or copyable by default, thats only nesseccary if you want to change the whole stored
 * object by assigning. If not, the data object can be uncopyable and it should be used by
 * retrieving it's reference with get-methods.
 *
 * If the property type is a standart c++ type like int or bool, the defualt value can't be set by using its
 * constructor. Therefore a interface for setting default values is added to the property. If you want
 * to assign a default value you just need to add a struct default_value which returns the wanted default
 * value with the operator(). If you don't want a default value, just don't add the struct. The implementation
 * assignes the default value to the property, therefore it should only be used with assignalble types.
 *
 *
 * A property implementation for storing integers at a graph edge with the identifier
 * 'test'property' may look like that:
 * @code
 * struct test_property {
 *  typedef int type;
 * }
 * @endcode
 *
 * The same property with a default value would be implemented like this:
 * @code
 * struct test_property {
 *  typedef int type;
 *  struct default_value {
 *           int operator()() {
 *               return 3;
 *           };
 *      };
 * }
 * @endcode
 * 
 * An additional usefull feature is the change tracking. It can be activated for every property individual by
 * supplying a \a change_tracking type in the property definition (either by typedef or by a new struct).
 * It is then possible to query if the property was changed (and reset the chang flag) and furthermore one 
 * can connect to a signal which is emited at property change. Our int property with change tracking enabled
 * would look like this:
 * 
 * @code
 * struct test_property {
 *  typedef int type;
 *  struct change_tracking {};
 * }
 * @endcode
 * 
 * If you want to use properties in your class you should derive from PropertyOwner class, as it doas all the
 * hanling needed and gives you get and set functions which work with the designed identifiers.
 * 
 * @code
 * ProeprtyOwner<mpl::vector<test_property>> owner;
 * int t = owner->get<test_property>();
 * owner->setProperty<test_property>(5);
 * @endcode
 * 
 *
 * @{ */

/**
 * @brief Exeption for property errors
 *
 * This exception is thrown when a property related error is detected, for example if a objects is ask for a
 * property which it does not own. This exceptions own the error-code range from 300-399.
 **/
struct property_error : virtual boost::exception { };

/**
 * @brief Signal for change of a property
 *
 * This signal is emmited when the property given as template parameter is changed.
 **/
template<typename Property>
struct onChange {};

namespace details {

/** @addtogroup Metafunctions
 * @{
 * @brief Type traits to detect if the property has a default value
 *
 * If the user want to provide a default value for a property than he adds a default_value struct with a
 * operator() which returns the default value. To check if the this struct is available we add a type traits
 * which searches for this special struct.
 */
BOOST_MPL_HAS_XXX_TRAIT_DEF(default_value)

/**
 * @brief Type traits to detect if the property shall be under change control
 *
 * If the user want to detect property changes and connect signals on such an event he need to specify it by
 * adding a typedef or a struct called change_tracking to the property. To check if the such a type is
 * available we add a type traits which searches for this special type.
 */
BOOST_MPL_HAS_XXX_TRAIT_DEF(change_tracking)


/**
 * @brief Metafunction to expose the property storage type
 **/
template<typename T>
struct property_type {
    typedef typename T::type type;
};

template<typename Prop, typename PropertyList>
struct has_property {
    typedef typename mpl::find<PropertyList, Prop>::type iterator;
    typedef typename mpl::not_<boost::is_same<iterator, typename mpl::end<PropertyList>::type>>::type type;
};

/**
 * @brief Metafunction to ensures that a property is in a lists
 *
 * Checks for the existence of a given Property in the mpl::vector supplied and adds the property if it is not found.
 * @tparam List mpl::vector with property types
 * @tparam Prop the property which should be in the supplied list
 **/
template<typename List, typename Prop>
struct ensure_property {

    typedef typename mpl::if_ <
    boost::is_same <
    typename mpl::find<List, Prop>::type,
             typename mpl::end<List>::type > ,
             typename mpl::push_back<List, Prop>::type,
             List >::type type;
};

/**
 * @brief Metafunction to ensures that multiple properties are in the lists
 *
 * Checks for the existence of all Property's given in the mpl::vector supplied and adds the properties if it is not found.
 * @tparam List mpl::vector with property types
 * @tparam PropList mpl::vector with properties which should be in the supplied list
 **/
template<typename List, typename PropList>
struct ensure_properties {

    typedef typename mpl::fold < PropList, List,
            mpl::if_ <
            boost::is_same< mpl::find<mpl::_1, mpl::_2>, mpl::end<mpl::_1> >,
            mpl::push_back<mpl::_1, mpl::_2>,
            mpl::_1
            > >::type type;
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

/**
 * @brief Property vector to a fusion sequence of bools
 *
 * If we want to track changes to a property we need a bool value to store the current state. Basicly we
 * only need the state for propertys with "chenge_control", but to keep the order for easy accessing we
 * genereate one state for every proeprty
 **/
template<typename T>
struct bs { //bool sequence
    template<typename T2>
    struct state_type {
        typedef bool type;
    };
    typedef typename mpl::transform<T, state_type<mpl::_1> >::type bv;
    typedef typename fusion::result_of::as_vector< bv >::type type;
};

/**
 * @brief Property vector to a fusion signal map of change events
 *
 * It creates a signal map which can be uses as input for SignalOwner class. Each proeprty gets a
 * onChange<Propert> signal
 **/
template<typename T>
struct sm { //signal map
    typedef typename mpl::fold<T, mpl::map<>, mpl::insert<mpl::_1, mpl::pair<onChange<mpl::_2>,
            boost::function1<void, property_type<mpl::_2> > > > >::type type;
};


/**@}*/

/**
 * @brief Functor to assign default values to property
 *
 * This functor holds a pointer to the PropertyOwner in question. The operator() get the properties which
 * hold a default value and assigns this value to the property the owner holds.
 */
template<typename PropertyOwner>
struct apply_default {

    PropertyOwner* owner;
    apply_default(PropertyOwner* o) : owner(o) {};
    template<typename T>
    void operator()(const T& t) {
        owner->template setProperty<T>(typename T::default_value()());
    };
};

/**
 * @brief Class for property handling in dcm
 *
 * To ease the work with properties this class is provided. It receives all the properties, which shall be
 * handled, in a mpl::vector typelist as its template argument. Than easy access to all properties by get
 * and set functions is achieved. Furthermore change tracking is available for individual properties and the
 * whole set. On changes a signal is emmited to whcih one can connect arbitrary functions.
 *
 **/
template<typename PropertyList>
struct PropertyOwner : public SignalOwner<typename details::sm<PropertyList>::type> {

    /**
    * @brief Constructor assigning default values
    *
    * It's important to initialise the property fusion sequences with this constructor
    * as much handling has to be done to ensure the users default values are added correctly
    **/
    PropertyOwner();

    /**
    * @brief Access properties
    *
    * Returns a reference to the propertys actual value. The property type has to be owned by this class,
    * which means it needs to be in the typelist that was given as template parameter to this class.
    * @tparam Prop property type which should be accessed
    * @return const Prop::type& a reference to the properties actual value.
    **/
    template<typename Prop>
    const typename Prop::type& getProperty() const;

    /**
       * @brief Set properties
       *
       * Sets the value of a specified property. The property type has to be owned by this class,
       * which means it needs to be in the typelist that was given as template parameter to this class. This
       * function is used
       * @tparam Prop property type which should be setProperty
       * @param value value of type Prop::type which should be set in this object
       **/
    template<typename Prop>
    typename boost::disable_if<details::has_change_tracking<Prop> >::type setProperty(const typename Prop::type& value);
    template<typename Prop>
    typename boost::enable_if<details::has_change_tracking<Prop> >::type setProperty(const typename Prop::type& value);

    /**
    * @brief Access properties non-const
    *
    * Don't use this unless abselutly nesseccary. It is impossible to detect if the reference was changed outside
    * of the owner, hence tracking needs to be updated manual. You should never ever store a refence to a property,
    * as changes can't be tracked either. This function is only available to comply with boost graph property maps 
    * and for properties whiche are to big to effieciently be copyed before and after change.
    * @remark Note that you can mark a property changed via this function with \ref markPropertyChanged
    * @tparam Prop property type which should be accessed
    * @return Prop::type& a reference to the properties actual value.
    **/
    template<typename Prop>
    typename Prop::type& getPropertyAccessible();

    //*********************************//
    // Functions for change management //
    //*********************************//

    /**
    * @brief Check if the property was changed after the last change acknowledgement
    *
    * If the property has been initiaised with a default value or any change was acknowledged it counts
    * as unchanged. However, if the property value was set before and not acknowledged it returns false
    * @tparam Prop property type which should checked for change
    * @return bool true if the property was changed after the last acknowledgement
    **/
    template<typename Prop>
    bool isPropertyChanged() const;

    /**
     * @brief Acknowledge property change
     *
     * Marks the property as unchanged. This can be used to notice that the change was processed.
     * @tparam Prop the property to mar as changed
     **/
    template<typename Prop>
    void acknowledgePropertyChange();

    /**
     * @brief Check if any property was changed
     * 
     * @return bool true if any property has the change flag set to true, i.e. isPropertyChanged() returns true
     */
    bool hasPropertyChanges() const;

    /**
     * @brief Acknowledge every property
     *
     * Sets the change flag for every property to false
     */
    void acknowledgePropertyChanges();
    
    /**
     * @brief Markes a property as changed 
     * Somesimte it is needed to mark a property as changed regardless of the tracked state. This for example happens
     * when one changes a property in an untrackable way via \ref getPropertyAccessible . Then the change state can be
     * manually fixed via this function.
     * 
     * @tparam Prop the property to mar as changed
     */
    template<typename Prop>
    void markPropertyChanged();


private:
    /* It's imortant to not store the properties but their types. These types are
     * stored and accessed as fusion vector.
     * */
    typedef typename details::pts<PropertyList>::type Properties;

    /* To track changes to properties we store boolean state variables
     * */
    typedef typename details::bs<PropertyList>::type States;

    Properties  m_properties;
    States      m_states;
};

template<typename PropertyList>
PropertyOwner<PropertyList>::PropertyOwner() {

    //get a vies of all types which have a default value
    typedef typename mpl::filter_view<PropertyList, details::has_default_value<mpl::_> >::type view;
    //set the default value
    details::apply_default<PropertyOwner> func(this);
    mpl::for_each<view>(func);

    //set all change states to false initialy
    fusion::for_each(m_states, boost::phoenix::arg_names::arg1 = false);

#if defined(BOOST_MPL_CFG_NO_HAS_XXX)
    throw property_error() <<  boost::errinfo_errno(1) << error_message("no default values supported");
#endif
};

/**
 * @brief Convienience spezialisation to ease interaction with system class
 *
 * If no property is supplied for a PropertyOwner derived class, a mpl::void_ type will be retrieved. To
 * remove the burdon of checking for that type in the class definition this spezialisation is supplied.
 **/
template<>
struct PropertyOwner<mpl::void_> {
    template<typename Prop>
    typename Prop::type& getProperty() {
        throw property_error() <<  boost::errinfo_errno(2) << error_message("unknown property type");
    };

    template<typename Prop>
    void setProperty(typename Prop::type value) {
        throw property_error() <<  boost::errinfo_errno(2) << error_message("unknown property type");
    };
};

template<typename PropertyList>
template<typename Prop>
const typename Prop::type& PropertyOwner<PropertyList>::getProperty() const {

    typedef typename mpl::find<PropertyList, Prop>::type iterator;
    typedef typename mpl::distance<typename mpl::begin<PropertyList>::type, iterator>::type distance;
    BOOST_MPL_ASSERT((mpl::not_<boost::is_same<iterator, typename mpl::end<PropertyList>::type > >));
    return fusion::at<distance> (m_properties);
};

template<typename PropertyList>
template<typename Prop>
typename Prop::type& PropertyOwner<PropertyList>::getPropertyAccessible() {

    typedef typename mpl::find<PropertyList, Prop>::type iterator;
    typedef typename mpl::distance<typename mpl::begin<PropertyList>::type, iterator>::type distance;
    BOOST_MPL_ASSERT((mpl::not_<boost::is_same<iterator, typename mpl::end<PropertyList>::type > >));
    return fusion::at<distance> (m_properties);
};

template<typename PropertyList>
template<typename Prop>
typename boost::disable_if<details::has_change_tracking<Prop> >::type PropertyOwner<PropertyList>::setProperty(const typename Prop::type& value) {

    typedef typename mpl::find<PropertyList, Prop>::type iterator;
    typedef typename mpl::distance<typename mpl::begin<PropertyList>::type, iterator>::type distance;
    BOOST_MPL_ASSERT((mpl::not_<boost::is_same<iterator, typename mpl::end<PropertyList>::type > >));
    fusion::at<distance> (m_properties) = value;
};

template<typename PropertyList>
template<typename Prop>
typename boost::enable_if<details::has_change_tracking<Prop> >::type PropertyOwner<PropertyList>::setProperty(const typename Prop::type& value) {

    typedef typename mpl::find<PropertyList, Prop>::type iterator;
    typedef typename mpl::distance<typename mpl::begin<PropertyList>::type, iterator>::type distance;
    BOOST_MPL_ASSERT((mpl::not_<boost::is_same<iterator, typename mpl::end<PropertyList>::type > >));
    fusion::at<distance>(m_properties) = value;
    //keep track of the changes
    fusion::at<distance>(m_states) = true;
    //emit signal to notify of the change
    SignalOwner<typename details::sm<PropertyList>::type>::template emitSignal< onChange<Prop> >(value);
};

template<typename PropertyList>
template<typename Prop>
bool PropertyOwner<PropertyList>::isPropertyChanged() const {

    typedef typename mpl::find<PropertyList, Prop>::type iterator;
    typedef typename mpl::distance<typename mpl::begin<PropertyList>::type, iterator>::type distance;
    BOOST_MPL_ASSERT((mpl::not_<boost::is_same<iterator, typename mpl::end<PropertyList>::type > >));
    return fusion::at<distance>(m_states);
};

template<typename PropertyList>
template<typename Prop>
void PropertyOwner<PropertyList>::acknowledgePropertyChange() {

    typedef typename mpl::find<PropertyList, Prop>::type iterator;
    typedef typename mpl::distance<typename mpl::begin<PropertyList>::type, iterator>::type distance;
    BOOST_MPL_ASSERT((mpl::not_<boost::is_same<iterator, typename mpl::end<PropertyList>::type > >));
    fusion::at<distance>(m_states) = false;
};

template<typename PropertyList>
bool PropertyOwner<PropertyList>::hasPropertyChanges() const {

    bool res = false;
    fusion::for_each(m_states, boost::phoenix::ref(res) = boost::phoenix::ref(res) || boost::phoenix::arg_names::arg1);
    return res;
};

template<typename PropertyList>
void PropertyOwner<PropertyList>::acknowledgePropertyChanges() {

    fusion::for_each(m_states, boost::phoenix::arg_names::arg1 = false);
};

template<typename PropertyList>
template<typename Prop>
void PropertyOwner<PropertyList>::markPropertyChanged() {

    typedef typename mpl::find<PropertyList, Prop>::type iterator;
    typedef typename mpl::distance<typename mpl::begin<PropertyList>::type, iterator>::type distance;
    BOOST_MPL_ASSERT((mpl::not_<boost::is_same<iterator, typename mpl::end<PropertyList>::type > >));
    fusion::at<distance>(m_states) = true;
};

//now create some standart properties
//***********************************

/**
 * @brief Dummy property
 **/
struct empty_prop {
    typedef int type;
};

/**@}*/ //Property
/**@}*/ //Core

}//details
}//dcm

#endif //DCM_PROPERTY_H
