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

#ifndef DCM_GEOMETRY_3D_H
#define DCM_GEOMETRY_3D_H

#include <opendcm/core/geometry.hpp>
#include <test/module3d/m3d.hpp>

#include <boost/fusion/include/at_c.hpp>

namespace fusion = boost::fusion;

namespace dcm {
namespace geometry {
 
template<typename Kernel, bool MappedType = true>
struct Direction3 : public Geometry<Kernel, 3, MappedType> {
    
    using Geometry::VectorD;
    fusion::vector<VectorD> ValueSequence = fusion::make_vector(create<VectorD>());
    
    VectorD& value() {
        return fusion::at_c<0>(ValueSequence);
    };
};
           

template<typename Kernel, mpl::bool_ MappedType = mpl::false_>
struct Point3 : Geometry<Kernel, 3, MappedType> {
               
    using Geometry::Vector;
    fusion::vector<Vector> values;
        
    Vector& value() {
        return fusion::at_c<0>(values);
    };
};

template<typename Kernel, mpl::bool_ MappedType = mpl::false_>
struct Line3 : Geometry<Kernel, 3, MappedType> {
               
    using Geometry::Vector;
    fusion::vector<Vector, Vector> values;
        
    Vector& point() {
        return fusion::at_c<0>(values);
    };
    Vector& direction() {
        return fusion::at_c<1>(values);
    };
};

} //geometry

namespace numeric {
    
    template<typename Kernel>
    struct Point3 : public Geometry<Kernel, geometry::Point3> {};
    
    struct PointOnLine3 : public DependendGeometry<Kernel, geometry::Point3, geometry::Line3> {
        
        
    };
    
} //numeric

namespace modell {
  
  struct XYZ {
    /*Modell XYZ: 
     * 0 = X;
     * 1 = Y;
     * 2 = Z;
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
      a.finalize(t);
    };
  };
  
  struct XYZ2 {
    /*Modell XYZ2: two xyz parts after each other 
     * 0 = X;
     * 1 = Y;
     * 2 = Z;
     * 3 = X dir;
     * 4 = Y dir;
     * 5 = Z dir;
     */    
    template<typename Scalar, typename Accessor, typename Vector, typename Type>
    void extract(Type& t, Vector& v) {
      Accessor a;
      v(0) = a.template get<Scalar, 0>(t);
      v(1) = a.template get<Scalar, 1>(t);
      v(2) = a.template get<Scalar, 2>(t);
      v(3) = a.template get<Scalar, 3>(t);
      v(4) = a.template get<Scalar, 4>(t);
      v(5) = a.template get<Scalar, 5>(t);
    }
    
    template<typename Scalar, typename Accessor, typename Vector, typename Type>
    void inject(Type& t, Vector& v) {
      Accessor a;
      a.template set<Scalar, 0>(v(0), t);
      a.template set<Scalar, 1>(v(1), t);
      a.template set<Scalar, 2>(v(2), t);
      a.template set<Scalar, 3>(v(3), t);
      a.template set<Scalar, 4>(v(4), t);
      a.template set<Scalar, 5>(v(5), t);
      a.finalize(t);
    };
  };
  
  struct XYZ2P {
    /*Modell XYZ2P: two xyz parts after each other and one parameter
     * 0 = X;
     * 1 = Y;
     * 2 = Z;
     * 3 = X dir;
     * 4 = Y dir;
     * 5 = Z dir;
     * 6 = Parameter
     */    
    template<typename Scalar, typename Accessor, typename Vector, typename Type>
    void extract(Type& t, Vector& v) {
      Accessor a;
      v(0) = a.template get<Scalar, 0>(t);
      v(1) = a.template get<Scalar, 1>(t);
      v(2) = a.template get<Scalar, 2>(t);
      v(3) = a.template get<Scalar, 3>(t);
      v(4) = a.template get<Scalar, 4>(t);
      v(5) = a.template get<Scalar, 5>(t);
      v(6) = a.template get<Scalar, 6>(t);
    }
    
    template<typename Scalar, typename Accessor, typename Vector, typename Type>
    void inject(Type& t, Vector& v) {
      Accessor a;
      a.template set<Scalar, 0>(v(0), t);
      a.template set<Scalar, 1>(v(1), t);
      a.template set<Scalar, 2>(v(2), t);
      a.template set<Scalar, 3>(v(3), t);
      a.template set<Scalar, 4>(v(4), t);
      a.template set<Scalar, 5>(v(5), t);
      a.template set<Scalar, 6>(v(6), t);
      a.finalize(t);
    };
  };
  
}

//dummy accessor
struct dummy_accessor {

    template<typename Scalar, int ID, typename T>
    Scalar get(T& t) {
        return 1;
    };
    template<typename Scalar, int ID, typename T>
    void set(Scalar value, T& t) {
        //TODO: throw
    };
    template<typename T>
    void finalize(T& t) {};
};

//dummy geometry traits for boost blank, wil bever be used
template<>
struct geometry_traits<boost::blank> {
    typedef tag::direction3D tag;
    typedef modell::XYZ modell;
    typedef dummy_accessor accessor;
};

} //geometry



#endif //GCM_GEOMETRY_3D_H
