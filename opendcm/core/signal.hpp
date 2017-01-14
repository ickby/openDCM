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

#include <boost/enable_shared_from_this.hpp>

namespace mpl = boost::mpl;
namespace fusion = boost::fusion;

namespace dcm {

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
    typename boost::enable_if<mpl::has_key<SigMap, S>, Connection>::type
    connectSignal(typename mpl::at<SigMap, S>::type function);

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
    typename boost::enable_if<mpl::has_key<SigMap, S>, void>::type
    disconnectSignal(Connection c);

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
    
    template <typename S, typename ...Arg>
    typename boost::enable_if<mpl::has_key<SigMap, S>, void>::type  emitSignal(Arg&... args);
};

/***************************************************************************************************************
 *
 *                      IMPELEMNTATION
 * 
 * *************************************************************************************************************/

template<typename SigMap>
SignalOwner<SigMap>::SignalOwner() : m_signal_count(0) {};

template<typename SigMap>
template<typename S>
typename boost::enable_if<mpl::has_key<SigMap, S>, Connection>::type
SignalOwner<SigMap>::connectSignal(typename mpl::at<SigMap, S>::type function)
{
    typedef typename mpl::find<sig_name, S>::type iterator;
    typedef typename mpl::distance<typename mpl::begin<sig_name>::type, iterator>::type distance;
    typedef typename fusion::result_of::value_at<Signals, distance>::type map_type;
    map_type& map = fusion::at<distance>(m_signals);
    map[++m_signal_count] = function;
    return m_signal_count;
};

template<typename SigMap>
template<typename S>
typename boost::enable_if<mpl::has_key<SigMap, S>, void>::type
SignalOwner<SigMap>::disconnectSignal(Connection c)
{
    typedef typename mpl::find<sig_name, S>::type iterator;
    typedef typename mpl::distance<typename mpl::begin<sig_name>::type, iterator>::type distance;

    typedef typename fusion::result_of::value_at<Signals, distance>::type map_type;
    map_type& map = fusion::at<distance>(m_signals);
    map.erase(c);
};

template<typename SigMap>
template <typename S, typename ...Arg>
typename boost::enable_if<mpl::has_key<SigMap, S>, void>::type
SignalOwner<SigMap>::emitSignal(Arg&... args) { 
    
    typedef typename mpl::find<sig_name, S>::type iterator;
    typedef typename mpl::distance<typename mpl::begin<sig_name>::type, iterator>::type distance;
    typedef typename fusion::result_of::value_at<Signals, distance>::type map_type;
    map_type& map = fusion::at<distance>(m_signals);
    for (typename map_type::iterator it=map.begin(); it != map.end(); it++)
        (it->second)(args...);
};

}; //details
}; //dcm

#ifndef DCM_EXTERNAL_CORE
//#include "imp/signal_imp.hpp"
#endif

#endif //DCM_SIGNAL_H


