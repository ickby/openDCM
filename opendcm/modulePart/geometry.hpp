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

#ifndef GCM_GEOMETRY_PART_H
#define GCM_GEOMETRY_PART_H

#include <opendcm/core/geometry.hpp>

namespace dcm {
namespace tag {

struct part  {};

}

namespace modell {
  
  struct quaternion_wxyz_vec3 {
    /*Modell XYZ: 
     * 0 = w;
     * 1 = x;
     * 2 = y;
     * 3 = z;
     */    
    template<typename Scalar, typename Accessor, typename Vector1, typename Vector2, typename Type>
    void extract(Type& t, Vector1& v, Vector2& v2) {
      //Vector is a Quaternion here
      Accessor a;
      v.w() = a.template get<Scalar, 0>(t);
      v.x() = a.template get<Scalar, 1>(t);
      v.y() = a.template get<Scalar, 2>(t);
      v.z() = a.template get<Scalar, 3>(t);
      //Vector2 is a Eigen::Vector3
      v2(0) = a.template get<Scalar, 4>(t);
      v2(1) = a.template get<Scalar, 5>(t);
      v2(2) = a.template get<Scalar, 6>(t);
      
    }
    
    template<typename Scalar, typename Accessor, typename Vector1, typename Vector2, typename Type>
    void inject(Type& t, Vector1& v, Vector2& v2) {
      Accessor a;
      a.template set<Scalar, 0>(v.w(), t);
      a.template set<Scalar, 1>(v.x(), t);
      a.template set<Scalar, 2>(v.y(), t);
      a.template set<Scalar, 3>(v.z(), t);
      
      a.template set<Scalar, 4>(v2(0), t);
      a.template set<Scalar, 5>(v2(1), t);
      a.template set<Scalar, 6>(v2(2), t);
      
      a.finalize(t);
    };
  };
}

}

#endif //GCM_GEOMETRY_PART_H
