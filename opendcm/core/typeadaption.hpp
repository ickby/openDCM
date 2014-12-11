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

#ifndef DCM_TYPEADAPTION_H
#define DCM_TYPEADAPTION_H

#include <boost/variant.hpp>

namespace dcm {
  

template< typename T>
struct geometry_traits {
    BOOST_MPL_ASSERT_MSG(false, NO_GEOMETRY_TRAITS_SPECIFIED_FOR_TYPE, (T));
};

namespace accessor {
    
struct OrderdBracket {

    template<typename Scalar, int ID, typename T>
    Scalar get(T& t) {
        return t[ID];
    };
    template<typename Scalar, int ID, typename T>
    void set(Scalar value, T& t) {
        t[ID] = value;
    };
    template<typename T>
    void finalize(T& t) {};
};

struct OrderdRoundbracket {

    template<typename Scalar, int ID, typename T>
    Scalar get(T& t) {
        return t(ID);
    };
    template<typename Scalar, int ID,  typename T>
    void set(Scalar value, T& t) {
        t(ID) = value;
    };
    template<typename T>
    void finalize(T& t) {};
};

}
    
}//dcm

#endif //DCM_TYPEADAPTION_H







