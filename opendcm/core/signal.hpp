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

#ifndef DCM_SIGNAL_H
#define DCM_SIGNAL_H

#include <iostream>
#include <map>

#include <boost/mpl/at.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/mpl/map.hpp>
#include <boost/mpl/key_type.hpp>
#include <boost/mpl/value_type.hpp>

#include <boost/fusion/include/as_vector.hpp>

#include <boost/preprocessor.hpp>
#include <boost/preprocessor/repetition/repeat.hpp>
#include <boost/preprocessor/cat.hpp>
#include <boost/preprocessor/repetition/enum.hpp>
#include <boost/preprocessor/repetition/enum_params.hpp>
#include <boost/preprocessor/repetition/enum_trailing_params.hpp>
#include <boost/preprocessor/repetition/enum_binary_params.hpp>

#include <boost/enable_shared_from_this.hpp>

#include "property.hpp"


namespace mpl = boost::mpl;
namespace fusion = boost::fusion;

/* Preprocessor implementation of emit signal. As we need many overloads with diffrent number of
 * templated parameters we use boost preprocessor to do the hard repetive work. The definition and
 * implementation are definded first as they need to be known before usage
 * */
#define EMIT_ARGUMENTS(z, n, data) \
    BOOST_PP_CAT(data, n)

#define EMIT_SIGNAL_CALL_DEF(z, n, data) \
    template < \
    typename S  \
    BOOST_PP_ENUM_TRAILING_PARAMS(n, typename Arg) \
    > \
    typename boost::enable_if<mpl::has_key<SigMap, S>, void>::Type \
    emitSignal( \
                     BOOST_PP_ENUM_BINARY_PARAMS(n, Arg, const& arg) \
                   );

namespace dcm {

//few standart signal names
struct remove {};

typedef int Connection;

namespace details {

/** @brief Class to handle signal management
 * 
 * **/
template<typename SigMap>
struct SignalOwner {
  
    SignalOwner();

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
    Connection connectSignal(typename mpl::at<SigMap, S>::type function);

    /**
    * @brief Disconnects a slot for a specific signal.
    *
    * Disconnects a slot so that it dosn't get called at signal emittion. It's important to
    * disconnect the slot by the same boost:function it was connected with.
    *
    * @tparam S the signal type of interest
    * @param c connection with which the slot was initialy connected
    * @return void
    **/
    template<typename S>
    void disconnectSignal(Connection c);

    //with no vararg templates before c++11 we need preprocessor to create the overloads of emit signal we need
    BOOST_PP_REPEAT(5, EMIT_SIGNAL_CALL_DEF, ~)

protected:
    /*signal handling
     * extract all signal types to allow index search (inex search on signal functions would fail as same
     * signatures are supported for multiple signals). Create std::vectors to allow multiple slots per signal
     * and store these vectors in a fusion::vector for easy access.
     * */
    typedef typename mpl::fold < SigMap, mpl::vector<>,
            mpl::push_back<mpl::_1, mpl::key_type<SigMap, mpl::_2> > >::type sig_name;
    typedef typename mpl::fold < SigMap, mpl::vector<>,
            mpl::push_back<mpl::_1, mpl::value_type<SigMap, mpl::_2> > >::type sig_functions;
    typedef typename mpl::fold < sig_functions, mpl::vector<>,
            mpl::push_back<mpl::_1, std::map<int, mpl::_2> > >::type sig_vectors;
    typedef typename fusion::result_of::as_vector<sig_vectors>::type Signals;

    Signals m_signals;
    int  m_signal_count;
};

}; //details
}; //dcm

#ifndef DCM_EXTERNAL_CORE
#include "imp/signal_imp.hpp"
#endif

#endif //DCM_SIGNAL_H


