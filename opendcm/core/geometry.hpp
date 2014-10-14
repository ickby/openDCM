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


#ifndef DCM_GEOMETRY_H
#define DCM_GEOMETRY_H

#include <vector>

#include "kernel.hpp"
#include "transformation.hpp"

namespace mpl = boost::mpl;

namespace dcm {

//signal we use for recalculation
struct recalculated {};

//all supported geometry types for easy access and comparison
namespace geometry {

namespace weight {

enum types {
    parameter = 0,
    direction,
    point,
    line,
    segment,
    circle,
    arcOfCircle,
    ellipse,
    arcOfEllipse,
    plane,
    cylinder,
    sphere
};

}; //weight

//we need types as identifiers
struct Parameter        : mpl::int_<weight::parameter> {};
struct Direction        : mpl::int_<weight::direction> {};
struct Point            : mpl::int_<weight::point> {};
struct Line             : mpl::int_<weight::line> {};
struct Circle           : mpl::int_<weight::circle> {};
struct ArcOfCircle      : mpl::int_<weight::arcOfCircle> {};
struct Ellipse          : mpl::int_<weight::ellipse> {};
struct ArcOfEllipse     : mpl::int_<weight::arcOfEllipse> {};
struct Plane            : mpl::int_<weight::plane> {};
struct Cylinder         : mpl::int_<weight::cylinder> {};
struct Sphere           : mpl::int_<weight::sphere> {};
    
}//namespace geometry

namespace details {

namespace numeric {
    
    /**
     * @brief Base class for numeric handling of all geometry types
     * 
     * This class is the common base for all geometry types. It therefore has all methods which are 
     * needed to calculate the value dependend on the current parameters and also the derivative for 
     * each parameter.
     */
    
template< typename Kernel >
struct Geometry {
    
    typedef Kernel::Scalar                                         Scalar;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1>               Vector;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>  Matrix;
    
    struct Derivative {
    
        const Parameter<Kernel>& parameter();
        const Vector&            value();
        
    private:
        Vector              m_value;
        Parameter<Kernel>*  m_parameter;
        
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    };
    
    typedef std::vector< Derivative >::const_iterator DerivativeIterator;
    
    //creation functions
    virtual void setup();
    
    //access functions
    const Vector&       value();
    DerivativeIterator  derivativesBegin();
    DerivativeIterator  derivativesEnd();
    
protected:
    Eigen::Map< Vector >        m_value;
    std::vector< Parameter* >   m_parameters;
    std::vector< Derivative* >  m_derivatives;    
};

template< typename Kernel >
struct DependendGeometry : Geometry<Kernel> {
  
    using Geometry::Scalar;
    using Geometry::Vector;
    
    void setBaseGeometry(Geometry<Kernel>* g);
    
protected:
    Vector              m_value;
    Geometry<Kernel>*   m_base;
};

template< typename Kernel, int Dimension >
struct ClusterGeometry : DependendGeometry<Kernel> {
    
    using DependendGeometry::Scalar;
    using DependendGeometry::Vector;
    typedef Transform<Scalar, Dimension> Transformation;
    
    void transform(const Transformation& transform);
    
protected:
    Transformation m_cumulated;
};
    
} //numeric
} //details
} //dcm

#ifndef DCM_EXTERNAL_CORE
#include "imp/geometry_imp.hpp"
#endif
#endif // GCM_GEOMETRY_H
