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

#include <boost/mpl/bool.hpp>
#include <boost/mpl/transform.hpp>
#include <boost/mpl/apply_wrap.hpp>
#include <boost/fusion/include/vector.hpp>
#include <boost/fusion/include/make_vector.hpp>
#include <boost/fusion/include/as_vector.hpp>
#include <boost/fusion/include/mpl.hpp>
#include <boost/graph/graph_concepts.hpp>

#include "kernel.hpp"
#include "transformation.hpp"

namespace mpl    = boost::mpl;
namespace fusion = boost::fusion;

namespace dcm {

/**
 * @brief Geometric primitives handling
 * 
 * A geometric primitive is a geometry type like a point in 2D or 3D, a line or whatever. Primitive in this
 * context decribes the fact that such a type is only used to store and access the data needed to fully
 * define the geometry. No data conversion or any calculation can be done. There are however two requirements
 * defined for a geometric primitive: First the values need to be stored in a fusion vector named \a m_storage. 
 * This is important for automatic handling of primitives regardless of type. Second it must be possible
 * to initiate a geometric primitive either as storage type (realy holding the values) or as pointer or map 
 * type, referencing the data from somewhere else. Furtherore it must comply with the dcm \ref Kernel 
 * infrastructure.
 * 
 * As an example let's create a point2d geometric primitive. We get the \ref Kernel as always template argument.
 * Furthermore we have a bool template parameter which describes if the point shall be a storage or a refernce.
 * 
 * \code{.cpp}
 * template<typename Kernel, bool Map>
 * struct point2D {
 *      typedef Kernel::Scalar Scalar;
 *      boost::fusion::vector<Scalar, Scalar> m_storage;
 * 
 *      Scalar& x() {return boost::fusion::at_c<0>(m_storage);
 *      Scalar& y() {return boost::fusion::at_c<1>(m_storage);
 * }
 * 
 * template<typename Kernel>
 * struct point2D<Kernel, true {
 *      typedef Kernel::Scalar Scalar;
 *      boost::fusion::vector<Scalar*, Scalar*> m_storage;
 * 
 *      Scalar* x() {return boost::fusion::at_c<0>(m_storage);
 *      Scalar* y() {return boost::fusion::at_c<1>(m_storage);
 * }
 * \endcode
 * 
 * The first requirement has been followed by creating a fusion vector called m_stroage which holds our two 
 * data entries. The second requeirement has been fullfilled for specializing the point type if a map type is
 * required and then replacing the Scalar values with pointers to those.
 * 
 * \note The fusion vector must have the same size and order for map and storage type. Furthermore the enries
 * in the map type primitive must be constructible from pointers to the storage type. 
 *
 * As the example above needs quite some boilerplate for a primitve type helper classes are provided in this 
 * namespace which ease the creation of geometric primitives. 
 */
namespace geometry {

    /**
     * @brief Storage types from which a geometric primitive can be build 
     * 
     * This namespace contains types which abstract the difference between a map and a storage and the 
     * difference in construction method. This makes it easy to describe the most common types independent
     * from the geometric primtive requeirements. For example a linear algebra vector of size 3 can be given
     * as storage::Vector<3>. This type hold all nesseccary boilerplate to define and construct map and 
     * storage types.
     */
namespace storage  {

    /**
     * @brief A simple scalar parameter
     * 
     * If you need a storage of one scalar value this is the class to use. For the storage it creates the 
     * \ref Kernel defined scalar and for a map it creates a pointer to such a scalar. The pointer will be
     * nullptr initialized.
     */
    struct Parameter {
    
    template<typename Kernel, bool Map>
    struct create {
        typedef typename Kernel::Scalar type;
        static type build() {
            return type();
        };
    };
    template<typename Kernel>
    struct create<Kernel, true> {
        typedef typename Kernel::Scalar* type;
        static type build() {
            return type(nullptr);
        };
    };
};
    
/**
 * @brief A two dimensionla \a Eigen matrix
 * 
 * A stroage for a full matrix, either fixed size (when \a Row or \Coloum are gives as integers > 0) or 
 * dynamical sized (when \a Row or \Coloum are gives as Eigen::Dynamic). For storage it creates a Eigen::Matrix
 * and for a map a Eigen::Map to such a matrix. The map will be nullptr initialized on construction.
 * 
 * \tparam Row The matrix's row count, or Eigen::Dynamic for dynamicaly sized rows.
 * \tparam Coloum The matrix's coloum count, or Eigen::Dynamic for dynamicaly sized coloumns.
 */
template<int Row, int Coloum>
struct Matrix {

    template<typename Kernel, bool Map>
    struct create {
        typedef Eigen::Matrix<typename Kernel::Scalar, Row, Coloum> type;
        static type build() {
            return type();
        };
    };
    template<typename Kernel>
    struct create<Kernel, true> {
        typedef Eigen::Map<Eigen::Matrix<typename Kernel::Scalar, Row, Coloum>> type;
        static type build() {
            return type(nullptr);
        };
    };
};

/**
 * @brief A row \a Eigen vector 
 * 
 * A stroage for a row vector, either fixed size (when \a Size is gives as integer > 0) or 
 * dynamical sized (when \a Size is gives as Eigen::Dynamic). For storage it creates a Eigen::Matrix
 * and for a map a Eigen::Map to such a matrix, both with the coloum size = 1. The map will be 
 * nullptr initialized on construction.
 * 
 * \tparam Size The vectors row count, or Eigen::Dynamic for a dynamicaly sized vector.
 */
template<int Size>
using Vector = Matrix<Size, 1>;

} //storages

/**
 * @brief Helper base struct for creating primitive gepometry types
 * 
 * Primitive geometry has to fullfill two requirements, the storage based on fusion::vector and the possibility
 * to be a map or storage type. To reduce boilerplate this class is given which handles the fusion::vector 
 * creation automativ for storage and map types dependend on the template parameter \a Map. To tell the struct
 * which types shall be stored (or mapped) it is possible to pass an arbitrary number of storage types as 
 * template parameters. For example a line could be constructed using this class in the following way:
 * \code
 * template<typename Kernel, bool MappedType = true>
 * struct Line3 : public geometry::Geometry<Kernel, MappedType,
 *            geometry::storage::Vector<3>, geometry::storage::Vector<3>> {
 *
 *    using geometry::Geometry<Kernel, MappedType, 
 *        geometry::storage::Vector<3>, geometry::storage::Vector<3>>::m_storage;
 *    
 *    auto point() -> decltype(fusion::at_c<0>(m_storage)) {
 *        return fusion::at_c<0>(m_storage);
 *    };
 *    
 *    auto direction() -> decltype(fusion::at_c<1>(m_storage)) {
 *        return fusion::at_c<1>(m_storage);
 *    };
 * };
 * \endcode
 *
 * \tparam Kernel The math \ref Kernel in use
 * \tparam Map boolean which states if this geometric primitive shall be a storage or a map
 * \tparam StorageTypes Variadic sequence of storage types
 */
template<typename Kernel, bool Map, typename... StorageTypes>
struct Geometry {

protected:
    typedef mpl::vector<StorageTypes...>                                               StorageTypeSequence;
    typedef mpl::vector< typename StorageTypes::template create<Kernel, Map>::type...> StorageSequence;
    typedef typename fusion::result_of::as_vector<StorageSequence>::type               Storage;

    Storage m_storage = fusion::make_vector(StorageTypes::template create<Kernel, Map>::build()...);
};

}//geometry

namespace numeric {

/**
 * @brief Base class for numeric handling of geometry types
 *
 * This class is the common base for all geometry types when used in numerical calculations. It is used to
 * ease the handling of the derivatives, as well as dependend geometry and cluster geometry.
 *
 * To allow the one-to-one substitution of dependend geometry and normal geometry in equations it is
 * important that both derive from the \a Base geometric primitive which holds the values of the geometry, 
 * for example a Vector3 for a point or a point and a direction for a line.
 *
 * \tparam Kernel The math \ref Kernel in use
 * \tparam Base The geometric primitive on which the numeric geometry is based on as non-specialized template type
 */
template< typename Kernel, template<class, bool> class Base >
struct Geometry : public Base<Kernel, true> {

    typedef typename Kernel::Scalar  Scalar;

    typedef std::pair<Base<Kernel, false>, Parameter<Kernel>*>  Derivative;
    typedef typename std::vector<Derivative>::iterator          DerivativeIterator;

    //access functions
    DerivativeIterator  derivativesBegin();
    DerivativeIterator  derivativesEnd();

    //sometimes it is possible to optimize constraint derivative calculation when we are
    //a normal geometry (where parameret == value). To enable those optimisations we need
    //an indicator if this situation is given
    bool isIndependend();

protected:
    bool m_independent = true;
    std::vector< Parameter<Kernel>* >  m_parameters;
    std::vector< Derivative >  m_derivatives;
};

template< typename Kernel, template<class, bool> class Base, template<class, bool> class DBase >
struct DependendGeometry : public Geometry<Kernel, Base> {

    typedef Geometry<Kernel, Base> Inherited;
    DependendGeometry(Geometry<Kernel, DBase>* dg = nullptr) : Geometry<Kernel, Base>(), m_base(dg) {

        Inherited::m_independent = false;
        //automated mapping of the numeric geometry values to the dependend real vectors

    };

    void setBaseGeometry(Geometry<Kernel, DBase>* g);

protected:
    int                       m_parameterCount; //the amount of extra parameters this type needs
    Base<Kernel, false>       m_value;
    Geometry<Kernel, DBase>*  m_base;
};

template< typename Kernel, template<class, bool> class Derived,
          template<class, bool> class Base, int Dimension >
struct ClusterGeometry : DependendGeometry<Kernel, Derived, Base> {

    typedef typename Kernel::Scalar               Scalar;
    typedef details::Transform<Scalar, Dimension> Transformation;

    void transform(const Transformation& transform);

protected:
    Transformation m_cumulated;
};

} //numeric

namespace symbolic {


} //symbolic

} //dcm

#ifndef DCM_EXTERNAL_CORE
//#include "imp/geometry_imp.hpp"
#endif
#endif // GCM_GEOMETRY_H
