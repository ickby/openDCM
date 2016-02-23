/*
    openDCM, dimensional constraint manager
    Copyright (C) 2015  Stefan Troeger <stefantroeger@gmx.net>

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

#ifndef DCM_MAPDATATYPE_H
#define DCM_MAPDATATYPE_H

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <boost/graph/graph_concepts.hpp>

#include "transformation.hpp"
#include "logging.hpp"

namespace dcm {
namespace details {


/**
 * @brief Datatype which allows advanced mapping with Eigen3 types
 * 
 * 
 */
template<typename Kernel>
struct MapType {

    typedef typename Kernel::Scalar Scalar;
    
    MapType() : m_pointer(&m_storage) {}
    MapType(const Scalar& val) : m_storage(val), m_pointer(&m_storage) {}
    MapType(const MapType<Kernel>& val) : m_storage(*val.m_pointer), m_pointer(&m_storage) {}
    
    void operator=(const Scalar& val) {
        *m_pointer = val;
    }
    
    void operator=(const MapType<Kernel>& val) {
        *m_pointer = val.m_storage;
    }
    
    template<typename K> 
    friend std::ostream& operator<< (std::ostream &out, const MapType<K> &type);

    operator Scalar() const { return *m_pointer; }
    
    //allow to redirect the mapping
    void mapTo(double* pointer) {
        m_pointer = pointer;
    }
        
    //redirect the map to the internal storage and set the given value
    void setLocalValue(double val) {
        m_pointer = &m_storage;
        m_storage = val;
    }
    
private:
    Scalar* m_pointer;
    Scalar  m_storage;
};

}//details

}//dcm

template<typename Kernel>
std::ostream& operator<< (std::ostream &out, const dcm::details::MapType<Kernel> &type)
{
    // Since operator<< is a friend of the Point class, we can access
    // Point's members directly.
    out << *type->m_pointer;
    return out;
}

namespace Eigen {
template<typename Kernel> struct NumTraits<dcm::details::MapType<Kernel>>
 : NumTraits<typename Kernel::Scalar> // permits to get the epsilon, dummy_precision, lowest, highest functions
{
  typedef dcm::details::MapType<Kernel> Real;
  typedef dcm::details::MapType<Kernel> NonInteger;
  typedef dcm::details::MapType<Kernel> Nested;
  enum {
    IsComplex = 0,
    IsInteger = 0,
    IsSigned = 1,
    RequireInitialization = 0,
    ReadCost = 1,
    AddCost = 1,
    MulCost = 1
  };
};
}

#endif //DCM_MAPDATATYPE_H






