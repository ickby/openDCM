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

#ifndef DCM_GEOMETRY_PART_H
#define DCM_GEOMETRY_PART_H

#include <opendcm/core/geometry.hpp>

namespace dcm {
namespace tag {

struct part  {};

}

namespace modell {
  
  struct quaternion_wxyz {
    /*Modell XYZ: 
     * 0 = w;
     * 1 = x;
     * 2 = y;
     * 3 = z;
     */    
    template<typename Scalar, typename Accessor, typename Vector, typename Type>
    void extract(Type& t, Vector& v) {
      Accessor a;
      v(0) = a.template get<Scalar, 0>(t);
      v(1) = a.template get<Scalar, 1>(t);
      v(2) = a.template get<Scalar, 2>(t);
    }
    
    template<typename Scalar, typename Accessor, typename Vector, typename Type>
    void inject(Type& t, Vector& v) {
      Accessor a;
      a.template set<Scalar, 0>(v(0), t);
      a.template set<Scalar, 1>(v(1), t);
      a.template set<Scalar, 2>(v(2), t);
    };
  };
}

}

#endif //DCM_GEOMETRY_PART_H
