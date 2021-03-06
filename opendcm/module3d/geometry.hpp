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

#ifndef DCM_GEOMETRY_3D_H
#define DCM_GEOMETRY_3D_H

#include <opendcm/core/geometry.hpp>
#include <boost/fusion/include/at_c.hpp>

namespace fusion = boost::fusion;

namespace dcm {

//the geometry primitives we handle in the 3d module
namespace geometry {

template<typename Kernel>
struct Point3 : public Geometry<Kernel, numeric::Vector<Kernel, 3>> {

    typedef typename Kernel::Scalar Scalar;
    typedef Geometry<Kernel, numeric::Vector<Kernel, 3>> Base;
    using Base::m_storage;

    auto point()->decltype(fusion::at_c<0>(m_storage)) {
        return fusion::at_c<0>(m_storage);
    };
    
    template<typename Derived>
    Point3<Kernel>& transform(const details::TransformBase<Derived>& t) {
        point() = t*point();
        return *this;
    };
    
    template<typename Derived>
    Point3<Kernel>  transformed(const details::TransformBase<Derived>& t) {
        Point3<Kernel> copy(*this);
        copy.transform(t);
        return copy;
    };
    
    static int id() {
        static int ID = Base::generateId();
        return ID;
    };
};

template<typename Kernel>
struct Line3 : public Geometry<Kernel, numeric::Vector<Kernel, 3>, numeric::Vector<Kernel, 3>> {

    typedef typename Kernel::Scalar Scalar;
    typedef Geometry<Kernel, numeric::Vector<Kernel, 3>, numeric::Vector<Kernel, 3>> Base;
    using Base::m_storage;

    auto point()->decltype(fusion::at_c<0>(m_storage)) {
        return fusion::at_c<0>(m_storage);
    };

    auto direction()->decltype(fusion::at_c<1>(m_storage)) {
        return fusion::at_c<1>(m_storage);
    };
    
    template<typename Derived>
    Line3<Kernel>& transform(const details::TransformBase<Derived>& t) {
        t.transform(point());
        t.rotate(direction());
        return *this;
    };
    
    template<typename Derived>
    Line3<Kernel>  transformed(const details::TransformBase<Derived>& t) {
        Line3<Kernel> copy(*this);
        copy.transform(t);
        return copy;
    };
    
    static int id() {
        static int ID = Base::generateId();
        return ID;
    };
};

template<typename Kernel>
struct Plane : public Geometry<Kernel, numeric::Vector<Kernel, 3>, numeric::Vector<Kernel, 3>> {

    typedef typename Kernel::Scalar Scalar;
    typedef Geometry<Kernel, numeric::Vector<Kernel, 3>, numeric::Vector<Kernel, 3>> Base;
    using Base::m_storage;

    auto point()->decltype(fusion::at_c<0>(m_storage)) {
        return fusion::at_c<0>(m_storage);
    };

    auto direction()->decltype(fusion::at_c<1>(m_storage)) {
        return fusion::at_c<1>(m_storage);
    };
    
    template<typename Derived>
    Plane<Kernel>& transform(const details::TransformBase<Derived>& t) {
        t.transform(point());
        t.rotate(direction());
        return *this;
    };
    
    template<typename Derived>
    Plane<Kernel>  transformed(const details::TransformBase<Derived>& t) {
        Plane<Kernel> copy(*this);
        copy.transform(t);
        return copy;
    };
    
    static int id() {
        static int ID = Base::generateId();
        return ID;
    };
};

template<typename Kernel>
struct Sphere: public Geometry<Kernel, numeric::Vector<Kernel, 3>, typename Kernel::Scalar> {

    typedef typename Kernel::Scalar Scalar;
    typedef Geometry<Kernel, numeric::Vector<Kernel, 3>, numeric::Vector<Kernel, 3>, Scalar> Base;
    using Base::m_storage;

    auto point()->decltype(fusion::at_c<0>(m_storage)) {
        return fusion::at_c<0>(m_storage);
    };

    Scalar& radius() {
        return Base::rmPtr(fusion::at_c<1>(m_storage));
    };
    
    template<typename Derived>
    Sphere<Kernel>& transform(const details::TransformBase<Derived>& t) {
        t.transform(point());
        return *this;
    };
    
    template<typename Derived>
    Sphere<Kernel>  transformed(const details::TransformBase<Derived>& t) {
        Sphere<Kernel> copy(*this);
        copy.transform(t);
        return copy;
    };
    
    static int id() {
        static int ID = Base::generateId();
        return ID;
    };
};

template<typename Kernel>
struct Cylinder : public Geometry<Kernel, numeric::Vector<Kernel, 3>, numeric::Vector<Kernel, 3>, typename Kernel::Scalar> {

    typedef typename Kernel::Scalar Scalar;
    typedef Geometry<Kernel, numeric::Vector<Kernel, 3>, numeric::Vector<Kernel, 3>, Scalar> Base;
    using Base::m_storage;

    auto point()->decltype(fusion::at_c<0>(m_storage)) {
        return fusion::at_c<0>(m_storage);
    };

    auto direction()->decltype(fusion::at_c<1>(m_storage)) {
        return fusion::at_c<1>(m_storage);
    };

    Scalar& radius() {
        return Base::rmPtr(fusion::at_c<2>(m_storage));
    };
    
    template<typename Derived>
    Cylinder<Kernel>& transform(const details::TransformBase<Derived>& t) {
        t.transform(point());
        t.rotate(direction());
        return *this;
    };
    
    template<typename Derived>
    Cylinder<Kernel>  transformed(const details::TransformBase<Derived>& t) {
        Cylinder<Kernel> copy(*this);
        copy.transform(t);
        return copy;
    };
    
    static int id() {
        static int ID = Base::generateId();
        return ID;
    };
};

}//geometry

//the user-exposed geometry types for use in the geometry traits
typedef geometry::Point3<DummyKernel>     Point3;
typedef geometry::Line3<DummyKernel>      Line3;
typedef geometry::Plane<DummyKernel>      Plane;
typedef geometry::Cylinder<DummyKernel>   Cylinder;
typedef geometry::Sphere<DummyKernel>     Sphere;

namespace modell {
    
struct CartesianDirection {
    /*Modell XYZ:
     * 0 = X;
     * 1 = Y;
     * 2 = Z;
     */
    template<typename Scalar, typename Accessor, typename Primitive, typename Type>
    void extract(Type& t, Primitive& v) {
        Accessor a;
        v.direction()(0) = a.template get<Scalar, 0>(t);
        v.direction()(1) = a.template get<Scalar, 1>(t);
        v.direction()(2) = a.template get<Scalar, 2>(t);
    }

    template<typename Scalar, typename Accessor, typename Primitive, typename Type>
    void inject(Type& t, Primitive& v) {
        Accessor a;
        a.template set<Scalar, 0>(v.direction()(0), t);
        a.template set<Scalar, 1>(v.direction()(1), t);
        a.template set<Scalar, 2>(v.direction()(2), t);
        a.finalize(t);
    };
};
   
struct CartesianPoint {
    /*Modell XYZ:
     * 0 = X;
     * 1 = Y;
     * 2 = Z;
     */
    template<typename Scalar, typename Accessor, typename Primitive, typename Type>
    void extract(Type& t, Primitive& v) {
        Accessor a;
        v.point()(0) = a.template get<Scalar, 0>(t);
        v.point()(1) = a.template get<Scalar, 1>(t);
        v.point()(2) = a.template get<Scalar, 2>(t);
    }

    template<typename Scalar, typename Accessor, typename Primitive, typename Type>
    void inject(Type& t, Primitive& v) {
        Accessor a;
        a.template set<Scalar, 0>(v.point()(0), t);
        a.template set<Scalar, 1>(v.point()(1), t);
        a.template set<Scalar, 2>(v.point()(2), t);
        a.finalize(t);
    };
};

struct CartesianPointDirection {
    /*Modell XYZ2: two xyz parts after each other
     * 0 = X;
     * 1 = Y;
     * 2 = Z;
     * 3 = X dir;
     * 4 = Y dir;
     * 5 = Z dir;
     */
    template<typename Scalar, typename Accessor, typename Primitive, typename Type>
    void extract(Type& t, Primitive& v) {
        Accessor a;
        v.point()(0) = a.template get<Scalar, 0>(t);
        v.point()(1) = a.template get<Scalar, 1>(t);
        v.point()(2) = a.template get<Scalar, 2>(t);
        v.direction()(0) = a.template get<Scalar, 3>(t);
        v.direction()(1) = a.template get<Scalar, 4>(t);
        v.direction()(2) = a.template get<Scalar, 5>(t);
    }

    template<typename Scalar, typename Accessor, typename Primitive, typename Type>
    void inject(Type& t, Primitive& v) {
        Accessor a;
        a.template set<Scalar, 0>(v.point()(0), t);
        a.template set<Scalar, 1>(v.point()(1), t);
        a.template set<Scalar, 2>(v.point()(2), t);
        a.template set<Scalar, 3>(v.direction()(0), t);
        a.template set<Scalar, 4>(v.direction()(1), t);
        a.template set<Scalar, 5>(v.direction()(2), t);
        a.finalize(t);
    };
};

struct CartesianPointDirectionRadius {
    /*Modell XYZ2P: two xyz parts after each other and one parameter
     * 0 = X;
     * 1 = Y;
     * 2 = Z;
     * 3 = X dir;
     * 4 = Y dir;
     * 5 = Z dir;
     * 6 = Parameter
     */
    template<typename Scalar, typename Accessor, typename Primitive, typename Type>
    void extract(Type& t, Primitive& v) {
        Accessor a;
        v.point()(0) = a.template get<Scalar, 0>(t);
        v.point()(1) = a.template get<Scalar, 1>(t);
        v.point()(2) = a.template get<Scalar, 2>(t);
        v.direction()(0) = a.template get<Scalar, 3>(t);
        v.direction()(1) = a.template get<Scalar, 4>(t);
        v.direction()(2) = a.template get<Scalar, 5>(t);
        v.radius() = a.template get<Scalar, 6>(t);
    }

    template<typename Scalar, typename Accessor, typename Primitive, typename Type>
    void inject(Type& t, Primitive& v) {
        Accessor a;
        a.template set<Scalar, 0>(v.point()(0), t);
        a.template set<Scalar, 1>(v.point()(1), t);
        a.template set<Scalar, 2>(v.point()(2), t);
        a.template set<Scalar, 3>(v.direction()(0), t);
        a.template set<Scalar, 4>(v.direction()(1), t);
        a.template set<Scalar, 5>(v.direction()(2), t);
        a.template set<Scalar, 6>(v.radius(), t);
        a.finalize(t);
    };
};

}

namespace accessor {
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
};

//dummy geometry traits for boost blank, wil bever be used
template<>
struct geometry_traits<boost::blank> {
    typedef Point3 type;
    typedef modell::CartesianPoint modell;
    typedef accessor::dummy_accessor accessor;
};

}//dcm

#endif //DCM_GEOMETRY_3D_H
