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

#ifndef NS2_OBJECT_H
#define NS2_OBJECT_H

#include <iostream>

#include <boost/mpl/at.hpp>
#include <boost/mpl/if.hpp>
#include <boost/mpl/void.hpp>

#include <boost/fusion/include/as_vector.hpp>
#include <boost/fusion/include/mpl.hpp>
#include <boost/fusion/include/at.hpp>

#include "property.hpp"

namespace mpl = boost::mpl;
namespace fusion = boost::fusion;

namespace dcm {

template<typename Sys, typename Obj>
struct Object {

public:
    template<typename Prop>
    typename Prop::type& getProperty() {
        typedef typename mpl::find<Properties, Prop>::type iterator;
        typedef typename mpl::distance<typename mpl::begin<Properties>::type, iterator>::type distance;
        return fusion::at<distance>(m_properties);
    };

    template<typename Prop>
    void setProperty(Prop& value) {
        typedef typename mpl::find<Properties, Prop>::type iterator;
        typedef typename mpl::distance<typename mpl::begin<Properties>::type, iterator>::type distance;
        fusion::at<distance>(m_properties) = value;
    };

protected:

    //typedef 
    typedef typename mpl::at<typename Sys::obj_properties, Obj>::type Mapped;
    typedef typename mpl::if_< typename boost::is_same<Mapped, mpl::void_>::type, mpl::vector<>, Mapped>::type Sequence;
    typedef typename fusion::result_of::as_vector<Sequence>::type Properties;

    Sys m_system;
    Properties m_properties;
};

}

#endif //NS2_OBJECT_H
