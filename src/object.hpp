/*
    openDCM, dimensional constraint manager
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

#ifndef DCM_OBJECT_H
#define DCM_OBJECT_H

#include <iostream>

#include <boost/mpl/at.hpp>
#include <boost/mpl/if.hpp>
#include <boost/mpl/void.hpp>
#include <boost/mpl/vector.hpp>

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

template<typename T>
void pretty(T t) {
  std::cout<<__PRETTY_FUNCTION__<<std::endl;
};

template<typename Sys, typename Obj, typename Sig>
struct Object {

    Object(Sys& system) : m_system(system) {};

    template<typename Prop>
    typename Prop::type& getProperty() {
        typedef typename mpl::find<Sequence, Prop>::type iterator;
        typedef typename mpl::distance<typename mpl::begin<Sequence>::type, iterator>::type distance;
        return fusion::at<distance>(m_properties);
    };

    template<typename Prop>
    void setProperty(typename Prop::type value) {
        typedef typename mpl::find<Sequence, Prop>::type iterator;
        typedef typename mpl::distance<typename mpl::begin<Sequence>::type, iterator>::type distance;
        fusion::at<distance>(m_properties) = value;
    };

    template<typename S>
    void connectSignal( typename mpl::at<Sig, S>::type function ) {
        typedef typename mpl::find<sig_name, S>::type iterator;
        typedef typename mpl::distance<typename mpl::begin<sig_name>::type, iterator>::type distance;
        fusion::at<distance>(m_signals).push_back(function);
    };

    template<typename S>
    void disconnectSignal( typename mpl::at<Sig, S>::type function ) {
        typedef typename mpl::find<sig_name, S>::type iterator;
        typedef typename mpl::distance<typename mpl::begin<sig_name>::type, iterator>::type distance;

        typedef typename fusion::result_of::at<Signals, distance>::type result;
        result& vec = fusion::at<distance>(m_signals);
        vec.erase(std::remove(vec.begin(), vec.end(), function), vec.end());
    };

protected:

    //properties
    typedef typename mpl::at<typename Sys::object_properties, Obj>::type Mapped;
    typedef typename mpl::if_< boost::is_same<Mapped, mpl::void_ >, mpl::vector<>, Mapped>::type Sequence;
    typedef typename mpl::transform<Sequence, details::property_type<mpl::_1> >::type Typesequence;
    typedef typename fusion::result_of::as_vector<Typesequence>::type Properties;

    //signal handling
    typedef typename mpl::fold< Sig, mpl::vector<>,
    mpl::push_back<mpl::_1, details::map_key<mpl::_2> > >::type sig_name;
    typedef typename mpl::fold< Sig, mpl::vector<>,
    mpl::push_back<mpl::_1, details::map_val<mpl::_2> > >::type sig_functions;
    typedef typename mpl::fold< sig_functions, mpl::vector<>,
    mpl::push_back<mpl::_1, std::vector<mpl::_2> > >::type sig_vectors;
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
      typedef typename fusion::result_of::value_at<Signals, distance>::type result; \
      result& vec = fusion::at<distance>(m_signals); \
      for (typename result::iterator it=vec.begin(); it != vec.end(); it++) \
	(*it)(BOOST_PP_ENUM(n, EMIT_ARGUMENTS, arg)); \
  };

    BOOST_PP_REPEAT(5, EMIT_CALL, ~)



};

}

#endif //DCM_OBJECT_H


