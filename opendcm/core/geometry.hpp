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

#include <Eigen/Core>
#include <Eigen/Dense>

#include <boost/mpl/bool.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/mpl/range_c.hpp>
#include <boost/mpl/for_each.hpp>
#include <boost/fusion/include/vector.hpp>
#include <boost/fusion/include/make_vector.hpp>
#include <boost/fusion/include/as_vector.hpp>
#include <boost/fusion/include/mpl.hpp>
#include <boost/fusion/include/for_each.hpp>

#include "kernel.hpp"
#include "transformation.hpp"

#ifdef DCM_DEBUG
#include "defines.hpp"
#endif

namespace mpl    = boost::mpl;
namespace fusion = boost::fusion;

namespace dcm {

/**
 * @brief Geometric primitives handling
 *
 * A geometric primitive is a geometry type like a point in 2D or 3D, a line or whatever. Primitive in this
 * context decribes the fact that such a type is only used to store and access the data needed to fully
 * define the geometry. The only operation on the data defined here is tranforming. There exist a few
 * requirements for a geometric primitive: 
 * First the values need to be stored in a fusion vector named \a m_storage. This is important for automatic 
 * handling of primitives regardless of type. 
 * Second it must be possible to initiate a geometric primitive either as storage type (realy holding the values) or as pointer or map
 * type, referencing the data from somewhere else.
 * Third it must implement the interface for in place transformation returning a transformed copy. This interface 
 * consists of two function: 
 * \code
 * this_type& transform  (const Eigen::Transform<Scalar, Dim, Eigen::AffineCompact>&) //in place transformation, return reference
 * this_type  transformed(const Eigen::Transform<Scalar, Dim, Eigen::AffineCompact>&) //return transformed copy
 *\endcode
 * Finally a geometric primitive must be copy-constructable and assignable by the same type in map and storage form.
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
 * 
 *      point2d<Kernel, Map>& transform  (const Eigen::Transform<Scalar, 2, Eigen::AffineCompact>& t)...
 *      point2d<Kernel, Map>  transformed(const Eigen::Transform<Scalar, 2, Eigen::AffineCompact>& t)...
 * }
 *
 * template<typename Kernel>
 * struct point2D<Kernel, true {
 *      typedef Kernel::Scalar Scalar;
 *      boost::fusion::vector<Scalar*, Scalar*> m_storage;
 *
 *      Scalar* x() {return boost::fusion::at_c<0>(m_storage);
 *      Scalar* y() {return boost::fusion::at_c<1>(m_storage);
 * 
 *      point2d<Kernel, Map>& transform  (const details::Transform<Scalar,2>& t)...
 *      point2d<Kernel, Map>  transformed(const details::Transform<Scalar,2>& t)...
 * }
 * }
 * \endcode
 *
 * The first requirement has been followed by creating a fusion vector called m_stroage which holds our two
 * data entries. The second requeirement has been fullfilled for specializing the point type if a map type is
 * required and then replacing the Scalar values with pointers to those. Number three is also satisfied, both 
 * spezialisations implement the transform interface.
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
 * \note A parameter will never be transformed
 */
struct Parameter {

    template<typename Kernel, bool Map>
    struct create {
        typedef typename Kernel::Scalar type;
        static type build() {
            return type(0);
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
 * A stroage for a full matrix, either fixed size (when \a Row or \Column are gives as integers > 0) or
 * dynamical sized (when \a Row or \Column are gives as Eigen::Dynamic). For storage it creates a Eigen::Matrix
 * and for a map a Eigen::Map to such a matrix. The map will be nullptr initialized on construction.
 * 
 * \note Matrices can be transformed in multiple ways. A <3,3> matrix can be rotated
 *
 * \tparam Row The matrix's row count, or Eigen::Dynamic for dynamicaly sized rows.
 * \tparam Column The matrix's column count, or Eigen::Dynamic for dynamicaly sized columns.
 */
template<int Row, int Column>
struct Matrix {

    template<typename Kernel, bool Map>
    struct create {
        typedef Eigen::Matrix<typename Kernel::Scalar, Row, Column> type;
        static type build() {
            return type::Zero();
        };
    };
    template<typename Kernel>
    struct create<Kernel, true> {
        typedef Eigen::Map<Eigen::Matrix<typename Kernel::Scalar, Row, Column>> type;
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
 * and for a map a Eigen::Map to such a matrix, both with the column size = 1. The map will be
 * nullptr initialized on construction.
 *
 * \note A vector can be transformed if its \a Size is equal the systems dimension. Three possibilities exist: 
 * 0 
 * 
 * \tparam Size The vectors row count, or Eigen::Dynamic for a dynamicaly sized vector.
 */
template<int Size>
using Vector = Matrix<Size, 1>;

} //storages

/**
 * @brief Helper base struct for creating primitive gepometry types
 *
 * Primitive geometry has to fullfill four requirements: the storage based on fusion::vector, the possibility
 * to be a map or storage type, the transform interface and the copy-constructability and assignability 
 * independend of map or storage type. To reduce boilerplate this class is given which handles the fusion::vector
 * creation automaticly for storage and map types dependend on the template parameter \a Map. To tell the struct
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
 * \note The transform interface is not implemented by this class. Thats because it is impossible to determine
 * how each storage entry needs to be transformed. Also the space dimension is not known. Therefore it is 
 * the derived classes responibility to implement that interface. 
 *
 * \tparam Kernel The math \ref Kernel in use
 * \tparam Map boolean which states if this geometric primitive shall be a storage or a map
 * \tparam StorageTypes Variadic sequence of storage types
 */
template<typename Kernel, bool Map, typename... StorageTypes>
struct Geometry {

    typedef mpl::vector<StorageTypes...>                                               StorageTypeSequence;
    typedef mpl::vector< typename StorageTypes::template create<Kernel, Map>::type...> StorageSequence;
    typedef typename fusion::result_of::as_vector<StorageSequence>::type               Storage;
    
protected:
    Storage m_storage = fusion::make_vector(StorageTypes::template create<Kernel, Map>::build()...);

    //helper function for parameters: remove the pointer if the parameter is represented as such
    //and in every case return a reference
    template<typename T>
    T& rmPtr(T* t) {
        return *t;
    };
    template<typename T>
    T& rmPtr(T& t) {
        return t;
    };
};

/**
 * @brief Adaptor for easy primitive geometry handling
 * 
 * As a primitive geometry is mostly used without the template arguments it may be very cumbersome to use it. 
 * Especialy if a fully defined type is needed. Therefore this adaptor is provided. It is fully defined by a 
 * primitive geomerty even without the template arguments and provides template aliases to generate a primtive.
 * 
 * \param Base the primitive geometry which should be wrapped
 */
template<template<class, bool> class Base>
struct adaptor {
    template<typename Kernel, bool b> using primitive    = Base<Kernel, b>;
    template<typename Kernel>         using mapped  = Base<Kernel, true>;
    template<typename Kernel>         using storage = Base<Kernel, false>;
    
    typedef Base<numeric::DummyKernel, true> placeholder;
};

/**
 * @brief Extractor of geometry base from initialized type
 * 
 * As a primitive geometry is mostly used without the template arguments it may be very cumbersome to use it. 
 * Especialy if a fully defined type is needed. Therefore this extractor is provided. It allows to extract 
 * the primitive type from fully initialized ones. This makes it possible to store or handle all primitive
 * geometries as wanted and still access it's base whenever needed. 
 * 
 * \param Base the primitive geometry which should be wrapped
 */
template<typename T>
struct extractor {};

template<template<class, bool> class Base, typename Kernel, bool b>
struct extractor<Base<Kernel, b>> {
    
    template<typename K, bool b_> using primitive = Base<K,b_>;
    typedef Base<Kernel, true>  mapped;
    typedef Base<Kernel, false> storage;
};

}//geometry

namespace numeric {

template<typename T>
void pretty(T t) {
    std::cout << __PRETTY_FUNCTION__ << std::endl;
};

/**
 * @brief Base class for numeric handling of geometry types
 *
 * This class is the common base for all geometry types when used in numerical calculations. It provides an
 * interface from which all needed numerical data can be accessed. Relevant are the geometry values and the
 * derivatives dependend on the parameters this geometry depends on. This is enough for every case, no matter
 * if the geometry depends on another one, is in a cluster or if values depend on parameters. Therefore this
 * class is used in equations to do calculations undependend of the actual state of the geometry and only
 * dependend by the type of the provided primitive.
 *
  \note The numeric geometry types do not yet support variable size geometric primitives. Please ensure such
 * are not used as \a Base template parameter.
 *
 * \tparam Kernel The math \ref Kernel in use
 * \tparam Base The geometric primitive on which the numeric geometry is based on as non-specialized 
 *              template type
 */
template< typename Kernel, template<class, bool> class Base >
struct Geometry : public Base<Kernel, true> {

    typedef Base<Kernel, true>                         Inherited;
    typedef typename Kernel::Scalar                    Scalar;
    typedef SystemEntry<Kernel>                        Parameter;
    typedef typename std::vector<Parameter>::iterator  ParameterIterator;
    typedef std::pair<Base<Kernel, false>, Parameter>  Derivative;
    typedef typename std::vector<Derivative>::iterator DerivativeIterator;

    //mpl trickery to get a sequence counting from 0 to the size of stroage entries
    typedef mpl::range_c<int,0,
            mpl::size<typename Inherited::StorageTypeSequence>::value> StorageRange;
                
    Geometry() {
        
        fusion::for_each(Inherited::m_storage, Counter(m_parameterCount));
    };
    
    //allow access to parameters and derivatives
    std::vector< Parameter >&  parameters()   {
        return m_parameters;
    };
    std::vector< Derivative >& derivatives() {
        return m_derivatives;
    };
    
    //allow to get the number of parameters this geometry offers. This function can be called before
    //the geometry was initialized
    int parameterCount() {
        return m_parameterCount;
    };

    //sometimes it is possible to optimize constraint derivative calculation when we are
    //a normal geometry (where parameret == value). To enable those optimisations we need
    //an indicator if this situation is given
    bool isIndependend() {
        return m_independent;
    };

    //the assumptions amde here are only valid for pure geometry, so we give the derived
    //types the chance to override the initialisaion behaviour
    virtual void init(LinearSystem<Kernel>& sys) {
#ifdef DCM_DEBUG
        dcm_assert(!m_init);
        m_init = true;
#endif
        //mpl trickery to get a sequence counting from 0 to the size of stroage entries
        typedef mpl::range_c<int,0,
                mpl::size<typename Inherited::StorageTypeSequence>::value> StorageRange;
        //now iterate that sequence so we can access all storage elements with knowing the position
        //we are at (that is important to access the correct derivative storage position too)
        mpl::for_each<StorageRange>(Initializer<typename Inherited::Storage>(sys,
                                    Inherited::m_storage, m_parameters, m_derivatives));
    };

protected:
#ifdef DCM_DEBUG
    bool m_init = false;
#endif
    bool m_independent = true;
    int  m_parameterCount = 0;
    std::vector< Parameter >  m_parameters;
    std::vector< Derivative > m_derivatives;

    template<typename StorageType, bool InitDerivative = true>
    struct Initializer {

        LinearSystem<Kernel>&    m_system;
        std::vector<Parameter>&  m_entries;
        std::vector<Derivative>& m_derivatives;
        StorageType&             m_storage;

        Initializer(LinearSystem<Kernel>& s, StorageType& st, std::vector<Parameter>& vec,
                    std::vector<Derivative>& der) : m_system(s), m_entries(vec),
            m_derivatives(der), m_storage(st) {};

        template<typename T>
        void operator()(T t) {

            //get the parameter
            auto& t1 = fusion::at<T>(m_storage);
            std::vector<Parameter> v = m_system.mapParameter(t1);

            //create and set derivatives
            for(int i=0; i<v.size(); ++i) {
                m_derivatives.emplace_back(Derivative(Base<Kernel, false>(),  v[i]));
                
                if(InitDerivative) {
                    auto& t2 = fusion::at<T>(m_derivatives.back().first.m_storage);
                    setOne(t2, i);
                }
            };

            //set parameter
            std::move(v.begin(), v.end(), std::back_inserter(m_entries));
        }

        template<typename T>
        void setOne(Eigen::MatrixBase<T>& m, int n) {
            m(n) = 1;
        };
        void setOne(Scalar& s, int n)               {
            s = 1;
        };
    };
    
    struct Counter {

        int& count;
        Counter(int& c) : count(c) {
            count = 0;
        }; 

        template<typename T>
        void operator()(T& t) const {
            
            count += size(t);
        }
        
        template<typename T>
        void operator()(T*& t) const {
            
            count += size(*t);
        }

        template<typename T>
        int size(Eigen::MatrixBase<T>& m) const {
            return m.rows()*m.cols();
        };
        int size(Scalar& s) const {
            return 1;
        };
    };

public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

/**
 * @brief Base class for geometry which calculates its value from a set of parameters
 * 
 * Not for every geometry each value is automaticly a system parameter. Often you have the case that certain
 * values must be calculated from a few parameter. Then the number to variate in numeric solving is not 
 * anymore the value but the parameters. To ease the handling of such a case this class is given. It allows 
 * to specify a independent parameter storage and ensures that the values are not mapped into the linear system
 * but the parameters are. As the values still need to be a map (to enable polymophism with basic \ref Geometry)
 * a object local storge m_storage has been created which the geometry value maps are mapped to. So if you 
 * calculate the geometry values from the parameters you can either write to the normal maps or the local 
 * storage.
 * 
 * \tparam Kernel The math \ref Kernel in use
 * \tparam Base The geometric primitive this numeric geometry represents
 * \tparam ParameterStorageTypes Any number of storage types which describe the parameters
 */
template< typename Kernel, template<class, bool> class Base, typename... ParameterStorageTypes>
struct ParameterGeometry : public Geometry<Kernel, Base> {

    typedef typename Kernel::Scalar Scalar;
    typedef Geometry<Kernel, Base>  Inherited;
    typedef typename geometry::Geometry<Kernel, true, ParameterStorageTypes...>::Storage ParameterStorage;
    
    ParameterGeometry() {

        fusion::for_each(m_parameterStorage, typename Inherited::Counter(Inherited::m_parameterCount));
    };
    
    /**
     * @brief Initialization of the geometry for calculation
     * 
     * This function must be called before any access to the geometry is made. It ensures that all maps are
     * valid. If the geometry is accessed before calling init a segfault will happen. In debug mode an assert 
     * is called in this situation, but in release no warning will occure. 
     * The geometry can only be initialized with a \ref LinearSystem as the parameters of the geometry are maped
     * into this lienar system. It is therefore highly important that the given \ref LinearSystem is used for
     * all calculations involing this geometry.
     * 
     * @param sys LinearSystem the geometry is initalized with
     * @return void
     */
    virtual void init(LinearSystem<Kernel>& sys) {
#ifdef DCM_DEBUG
        dcm_assert(!Inherited::m_init);
        Inherited::m_init = true;
#endif
        //mpl trickery to get a sequence counting from 0 to the size of stroage entries
        typedef mpl::range_c<int,0,
                mpl::size<typename Inherited::StorageTypeSequence>::value> StorageRange;
        
        //now iterate that sequence so we can access all storage elements with knowing the position
        //we are at (that is important to access the correct derivative storage position too)
        mpl::for_each<StorageRange>(Remapper(Inherited::m_storage, m_value));
        
        //finally map the parameters and create derivatives for them (zero initialized)
        typedef mpl::range_c<int,0, mpl::size<ParameterStorage>::value> ParameterRange;
        mpl::for_each<ParameterRange>(typename Inherited::template Initializer<ParameterStorage, false>(
                                    sys, m_parameterStorage, Inherited::m_parameters, 
                                    Inherited::m_derivatives));
    };

protected:
    ParameterStorage     m_parameterStorage = fusion::make_vector(
                                              ParameterStorageTypes::template create<Kernel, true>::build()...);
    Base<Kernel, false>  m_value;

    struct Remapper {

        typename Inherited::Storage& m_map;
        Base<Kernel, false>&         m_storage;

        Remapper(typename Inherited::Storage& map, Base<Kernel, false>& storage)
            : m_map(map), m_storage(storage) {};

        template<typename T>
        void operator()(T t) {
            remap(fusion::at<T>(m_storage.m_storage), fusion::at<T>(m_map));
        };

        template<typename T>
        void remap(T& storage, Eigen::Map<T>& map) {
            new(&map) Eigen::Map<T>(&storage(0));
        };

        void remap(Scalar& storage, Scalar*& map) {
            map = &storage;
        };
    };
};

/**
 * @brief Base class for numeric geometry calculations of interdependend geomtries
 * 
 * Often a geometry can be calculated in relation to annother one. This reduces the amount of 
 * existing parameters. This base class eases the handling of such a case by allowing to specify
 * the base geometry and holding it for easy access. As it derives from ParameterGeometry the value
 * of \ref this can be calculated from a certain set of parameters and the other geometry only.
 * 
 * \tparam Kernel The math \ref Kenel in use
 * \tparam Base The geometric primitive this numeric geometry is based in
 * \tparam DBase The geometric primitive this numeric geometry depends on
 * \tparam ParameterStorageTypes Any number of storage types which describe the parameters
 */
template< typename Kernel, template<class, bool> class Base,
          template<class, bool> class DBase, typename... ParameterStorageTypes>
struct DependendGeometry : public ParameterGeometry<Kernel, Base, ParameterStorageTypes...>  {
    
    typedef ParameterGeometry<Kernel, Base, ParameterStorageTypes...> Inherited;
    typedef typename Geometry<Kernel, DBase>::Derivative              DependendDerivative;
    
    /**
     * @brief Setup the geometry \a this depends on
     * 
     * Allows to set the base geometry which is needed to calculate \a this geometrys value. It will 
     * automaticly setup new derivatives dependend on the base derivatives, however, the values are
     * not quranteed to have any value. They therefore need to be overridden before use.
     * \note This function should only be called if \a this and \a base have been initialized.
     * 
     * \param base The base geometry \a this depends on
     * @return void
     */
    void setBaseGeometry(Geometry<Kernel, DBase>* base) {
        
        dcm_assert(Inherited::m_init);
        
        m_base = base;
        //we are a dependend geometry, therefore our value depends on all parameters of base. 
        Inherited::m_parameters.insert(Inherited::m_parameters.end(),
                                       base->parameters().begin(), base->parameters().end());
        
        //This also means that we have a derivative for every parameter of base too.
        for(const typename Geometry<Kernel, DBase>::Derivative& d : base->derivatives())
                Inherited::m_derivatives.push_back(std::make_pair(Base<Kernel, false>(), d.second));
    };
    
protected:
    Geometry<Kernel, DBase>* m_base = nullptr;
};

} //numeric

namespace symbolic {

struct Geometry {
    
    int type;
};

template<typename Kernel, template<class, bool> class G>
struct TypeGeometry : public Geometry {

    typedef G<Kernel, false> PrimitiveGeometry; 
    
    PrimitiveGeometry& getPrimitveGeometry() {return m_geometry;};
    
    template<bool mapped>
    void setPrimitiveGeometry(const G<Kernel, mapped>& g) {
        m_geometry = g;
    };    
    
    void setGeometryID(int id) {
        type = id;
    }
    
protected:   
    PrimitiveGeometry m_geometry;
};

struct GeometryProperty {
    typedef Geometry* type;
    struct default_value {
        Geometry* operator()() {
            return nullptr;
        };
    };
    struct change_tracking{};
};

} //symbolic
} //dcm

#endif // DCM_GEOMETRY_H
