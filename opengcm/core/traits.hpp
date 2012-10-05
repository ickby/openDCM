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

#ifndef GCM_TRAITS_H
#define GCM_TRAITS_H

#include <boost/type_traits.hpp>
#include <boost/mpl/if.hpp>
#include <boost/mpl/void.hpp>

#include <string.h>

namespace mpl = boost::mpl;

namespace gcm {

template< typename T >
struct system_traits {
    typedef typename T::Kernel  Kernel;
    typedef typename T::Cluster Cluster;

    template<typename M>
    struct getModule {

        typedef typename mpl::if_< boost::is_base_of<M, typename T::Type1>, typename T::Type1, typename T::Type2 >::type test1;
        typedef typename mpl::if_< boost::is_base_of<M, test1>, test1, typename T::Type3 >::type test2;
        typedef typename mpl::if_< boost::is_base_of<M, test1>, test2, mpl::void_ >::type type;

        typedef boost::is_same<type, mpl::void_> has_module;
    };
};

template<typename T>
struct compare_traits {

    BOOST_MPL_ASSERT_MSG((boost::is_same<T, const char*>::value),
                         YOU_SHOULD_NOT_USE_THIS_TYPE_AS_IDENTIFIER,
                         (const char*));

    static bool compare(T& first, T& second) {
        return first == second;
    };
};

template<>
struct compare_traits<std::string> {

    static bool compare(std::string& first, std::string& second) {
        return !(first.compare(second));
    };
};

}

#endif //GCM_TRAITS_H