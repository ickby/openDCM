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

#ifndef GCM_GEOMETRY_PART_H
#define GCM_GEOMETRY_PART_H

#include <opendcm/core/geometry.hpp>
#include <opendcm/core/kernel.hpp>
#include <opendcm/module3d/cluster.hpp>

namespace dcm {

namespace modell {
  
  struct quaternion_wxyz_vec3 {
    /*Modell XYZ: 
     * 0 = w;
     * 1 = x;
     * 2 = y;
     * 3 = z;
     */    
    template<typename Scalar, typename Accessor, typename Primitive, typename Type>
    void extract(Type& t, Primitive& v) {
      
        Accessor a;
        dcm::details::Transform<Scalar, 3> trans;
        
        typename dcm::details::Transform<Scalar, 3>::Rotation Q( a.template get<Scalar, 0>(t),
                                                                 a.template get<Scalar, 1>(t),
                                                                 a.template get<Scalar, 2>(t),
                                                                 a.template get<Scalar, 3>(t));
        
        typename dcm::details::Transform<Scalar, 3>::Translation T( a.template get<Scalar, 4>(t),
                                                                    a.template get<Scalar, 5>(t),
                                                                    a.template get<Scalar, 6>(t));
        
        v.transform() =  dcm::details::Transform<Scalar, 3>(Q,T);
    }
    
    template<typename Scalar, typename Accessor, typename Primitive, typename Type>
    void inject(Type& t, Primitive& v) {
      
      typedef typename Eigen::Quaternion<Scalar>   Rotation;
      typedef typename Eigen::Matrix<Scalar, 3,1>  Translation;
      
      Accessor a;
      
      const Rotation& r(v.rotation());
      a.template set<Scalar, 0>(r.w(), t);
      a.template set<Scalar, 1>(r.x(), t);
      a.template set<Scalar, 2>(r.y(), t);
      a.template set<Scalar, 3>(r.z(), t);
      
      const Translation& tr(v.translation());
      a.template set<Scalar, 4>(tr(0), t);
      a.template set<Scalar, 5>(tr(1), t);
      a.template set<Scalar, 6>(tr(2), t);
      
      a.finalize(t);
    };
  };
}

//the geometry primitives we handle in the part odule
namespace geometry {

//A Part is just a fancy name for a cluster
template<typename Kernel>
struct Part3 : public dcm::geometry::Cluster3<Kernel> {};

}//geometry

//the user-exposed geometry types for use in the geometry traits
typedef dcm::geometry::Part3<DummyKernel> Part3;

}//dcm

#endif //GCM_GEOMETRY_PART_H
