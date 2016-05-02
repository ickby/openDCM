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
#include "equations.hpp"
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
 * Second it must implement the interface for in place transformation and returning a transformed copy. T
 * his interface consists of two function: 
 * \code
 * //in place transformation, return reference
 * this_type& transform  (const Eigen::Transform<Scalar, Dim, Eigen::AffineCompact>&)
 * //return transformed copy
 * this_type  transformed(const Eigen::Transform<Scalar, Dim, Eigen::AffineCompact>&) 
 *\endcode
 * Third and finale requirement is that a geometric primitive must be default-constructable, 
 * copy-constructable and assignable.
 * 
 * As an example let's create a point2d geometric primitive. We get the \ref Kernel as template argument. *
 * \code{.cpp}
 * template<typename Kernel>
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
 * \endcode
 *
 * The first requirement has been followed by creating a fusion vector called m_stroage which holds our two
 * data entries. The second requeirement is fullfilled by the fully provided transform interface. The third copy 
 * and assignable condition is meet by standart compiler-generated operators and constructors. 
 *
 * As the example above needs quite some boilerplate for a primitve type helper classes are provided in this
 * namespace which ease the creation of geometric primitives.
 */
namespace geometry {

/**
 * @brief Helper base struct for creating primitive gepometry types
 *
 * Primitive geometry has to fullfill three requirements: the storage based on fusion::vector, the 
 * transform interface and the copy-constructability and assignability. To reduce boilerplate this class
 * is given which handles the fusion::vector creation automaticly. To tell the struct
 * which types shall be stored it is possible to pass an arbitrary number of types as
 * template parameters. For example a line could be constructed using this class in the following way:
 * \code
 * template<typename Kernel>
 * struct Line3 : public geometry::Geometry<Kernel, Eigen::Vector3d, Eigen::Vector3d> {
 *
 *    typedef geometry::Geometry<Kernel, Eigen::Vector3d, Eigen::Vector3d> Base;
 *
 *    Eigen::Vector3d& point() {
 *        return fusion::at_c<0>(Base::m_storage);
 *    };
 *
 *    auto direction() -> decltype(fusion::at_c<1>(Base::m_storage)) {
 *        return fusion::at_c<1>(Base::m_storage);
 *    };
 * 
 *    Line3<Kernel>& transform  (const Eigen::Transform<Scalar, 3, Eigen::AffineCompact>& t)...
 *    Line3<Kernel>  transformed(const Eigen::Transform<Scalar, 3, Eigen::AffineCompact>& t)...        
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
template<typename Kernel, typename... StorageTypes>
struct Geometry {

    typedef mpl::vector< StorageTypes... >                               StorageSequence;
    typedef typename fusion::result_of::as_vector<StorageSequence>::type Storage;
    
protected:
    Storage m_storage;
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

template<template<class> class Base, typename Kernel>
struct extractor<Base<Kernel>> {
    
    template<typename K> using  primitive = Base<K>;
};

}//geometry

namespace detail {
//helper classes for numeric geometry
template<typename Kernel, typename StorageType, typename Equation, bool InitDerivative = true>
struct Initializer {

    typedef typename Kernel::Scalar Scalar;
    
    numeric::LinearSystem<Kernel>&                  m_system;
    std::vector<typename Equation::Parameter>&      m_entries;
    std::vector<typename Equation::DerivativePack>& m_derivatives;
    StorageType&                                    m_storage;

    Initializer(numeric::LinearSystem<Kernel>& s, StorageType& st, std::vector<typename Equation::Parameter>& vec,
                std::vector<typename Equation::DerivativePack>& der) : m_system(s), m_entries(vec),
        m_derivatives(der), m_storage(st) {};

    template<typename T>
    void operator()(T t) {

        //get the parameter
        auto& t1 = fusion::at<T>(m_storage);
        std::vector<typename Equation::Parameter> v = map(t1);

        //create and set derivatives
        for(int i=0; i<v.size(); ++i) {
            m_derivatives.emplace_back(typename Equation::DerivativePack(typename Equation::OutputType(),  v[i]));
            
            if(InitDerivative) {
                auto& t2 = fusion::at<T>(m_derivatives.back().first.m_storage);
                setOne(t2, i);
            }
        };

        //set parameters
        std::move(v.begin(), v.end(), std::back_inserter(m_entries));
    }

    template<typename T>
    void setOne(Eigen::MatrixBase<T>& m, int n) {
        m(n) = 1;
    };
    void setOne(Scalar& s, int n)               {
        s = 1;
    };
    
    template<typename T>
    std::vector<typename Equation::Parameter> map(const Eigen::MatrixBase<T>& m) const {
        std::vector<typename Equation::Parameter> vec(m.rows()*m.cols());
        
        for(int i=0; i<(m.rows()*m.cols()); ++i) 
            vec[i] = m_system.mapParameter();
        
        return vec;
    };
    
    std::vector<typename Equation::Parameter> map(const Scalar& s) const {
        std::vector<typename Equation::Parameter> vec(1);
        vec[0] = m_system.mapParameter();
        
        return vec;
    };
};

template<typename Kernel>
struct Counter {

    typedef typename Kernel::Scalar Scalar;
    
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

template<typename Equation>
struct Assigner {

    int& m_count;
    std::vector<typename Equation::Parameter>& m_params;

    Assigner(std::vector<typename Equation::Parameter>& p, int& c) : m_params(p), m_count(c) {
        m_count = -1;
    };

    template<typename T>
    void operator()(Eigen::MatrixBase<T>& t) const {
        for(int i=0; i<(t.rows()*t.cols()); ++i)
            t(i) = m_params[++m_count];
    }
    
    void operator()(typename Equation::KernelType::Scalar& t) const {
        t = m_params[++m_count];
    }
};

};//detail
    
namespace numeric {

template<typename Kernel, int i> 
using Vector = Eigen::Matrix<typename Kernel::Scalar, i, 1>;

template<typename Kernel, int i, int j> 
using Matrix = Eigen::Matrix<typename Kernel::Scalar, i, j>;

/**
 * @brief Base class for numeric handling of geometry types
 *
 * This class is used for numeric evaluation of primitive geometry. It does wrap the geometries and 
 * exposes them as equation for further use in complex mathematical systems. It is responsible for 
 * creating free parameters representing the primitive as well as initialising them correctly. Furthermore
 * it also exposes the correct derivatives for the free parameters. 
 * The intended use of this class is to give a simple one-line wrapper for all possible geometries where 
 * every number is used as a free parameter.
 *
  \note The numeric geometry types do not yet support variable size geometric primitives. Please ensure such
 * are not used as \a Base template parameter.
 *
 * \tparam Kernel The math \ref Kernel in use
 * \tparam Base The geometric primitive on which the numeric geometry is based on as non-specialized 
 *              template type
 */
template< typename Kernel, template<class> class Base >
struct Geometry : public Equation<Kernel, Base<Kernel>> {

    typedef Equation<Kernel, Base<Kernel>> Inherited;
    typedef typename Kernel::Scalar        Scalar;

    //mpl trickery to get a sequence counting from 0 to the size of stroage entries
    typedef mpl::range_c<int,0,
            mpl::size<typename Inherited::StorageSequence>::value> StorageRange;
                
    Geometry() {
        Inherited::m_complexity = Complexity::Complex;
        fusion::for_each(Inherited::m_storage, detail::Counter<Kernel>(Inherited::m_parameterCount));
    };
    

    //sometimes it is possible to optimize constraint derivative calculation when we are
    //a normal geometry (where parameret == value). To enable those optimisations we need
    //an indicator if this situation is given
    bool isIndependend() {
        return m_independent;
    };

    //the assumptions made here are only valid for pure geometry, so the derived
    //types must override the initialisaion behaviour
    virtual void init(LinearSystem<Kernel>& sys) {
#ifdef DCM_DEBUG
        dcm_assert(!Inherited::m_init);
        Inherited::m_init = true;
#endif
        //mpl trickery to get a sequence counting from 0 to the size of stroage entries
        typedef mpl::range_c<int,0,
                mpl::size<typename Inherited::StorageSequence>::value> StorageRange;
        //now iterate that sequence so we can access all storage elements with knowing the position
        //we are at (that is important to access the correct derivative storage position too)
        mpl::for_each<StorageRange>(detail::Initializer<Kernel, typename Inherited::Storage, Inherited>(sys,
                                    Inherited::m_storage, Inherited::m_parameters, 
                                    Inherited::m_derivatives));
    };
    
    //we actually do not really need to calculate anything, but we need to make sure the mapped 
    //values are move over to the output
    CALCULATE() {
        
        fusion::for_each(Inherited::m_storage, detail::Assigner<Inherited>(Inherited::m_parameters, m_counter));    
    };
    
protected:
    bool m_independent = true;
    
protected:
    int m_counter = -1; //we need a counter for every calculate, and we do not want a memory allocation 
                        //for every recalculate
};

/**
 * @brief Base class for geometry which calculates its value from a set of parameters
 * 
 * Not for every geometry each value is automaticly a system parameter. Often you have the case that certain
 * values must be calculated from a few parameters, and this class makes it easy to expose those cases 
 * as equations to the numerical solving process. The number to variate in numeric solving is not 
 * anymore the value of the output but the free parameters only. Hence this class allows to specify a 
 * independent parameter storage and ensures that the values of the result geometry are not 
 * mapped into the linear system but the parameters of the storage are. To specify the parameter 
 * storage the template parameter \tparam ParameterStorageTypes can be used.
 * 
 * \tparam Kernel The math \ref Kernel in use
 * \tparam Base The geometric primitive this numeric geometry represents
 * \tparam ParameterStorageTypes Any number of storage types which describe the parameters
 */
template< typename Kernel, template<class> class Base, typename... ParameterStorageTypes>
struct ParameterGeometry : public Equation<Kernel, Base<Kernel>> {

    typedef typename Kernel::Scalar         Scalar;
    typedef Equation<Kernel, Base<Kernel>>  Inherited;
    typedef typename geometry::Geometry<Kernel, ParameterStorageTypes...>::Storage ParameterStorage;
    
    ParameterGeometry() {
        fusion::for_each(m_parameterStorage, detail::Counter<Kernel>(Inherited::m_parameterCount));
    };
    
    //make sure the parameter storage is used, not the geometry one, for initialisation
    virtual void init(LinearSystem<Kernel>& sys) {
#ifdef DCM_DEBUG
        dcm_assert(!Inherited::m_init);
        Inherited::m_init = true;
#endif
        //mpl trickery to get a sequence counting from 0 to the size of stroage entries
        typedef mpl::range_c<int,0, mpl::size<ParameterStorage>::value> StorageRange;
        
        //now iterate that sequence so we can access all storage elements with knowing the position
        //we are at (that is important to access the correct derivative storage position too)
        mpl::for_each<StorageRange>(detail::Initializer<Kernel, ParameterStorage, Inherited>(sys,
                                    m_parameterStorage, Inherited::m_parameters, 
                                    Inherited::m_derivatives));
    };
    
    CALCULATE() {
        fusion::for_each(m_parameterStorage, 
                         detail::Assigner<Inherited>(Inherited::m_parameters, m_counter));    
    };

protected:
    ParameterStorage     m_parameterStorage;

    int m_counter = -1; //we need a counter for every calculate, and we do not want a memory allocation 
                        //for every recalculate
};

/**
 * @brief Base class for numeric geometry calculations of interdependend geomtries
 * 
 * Often a geometry can be calculated in relation to annother one. This reduces the amount of 
 * existing parameters. This base class eases the handling of such a case by deriving from UnaryEquation.
 * Hence it can hold a input equation, as all numeric geometrie wrappers are equations. This does of 
 * course also allow to calculate the input from any other equation, it must not be a geometry one. 
 * Equal to the \ref ParameterGeometry this class allows to specify a extra parameter storage which is 
 * used to define free parameters. Hence this equation type calculates its result from a input equation
 * and free parameters.
 * 
 * \tparam Kernel The math \ref Kenel in use
 * \tparam Input  The geometric primitive this numeric geometry is based in
 * \tparam Output The geometric primitive this numeric geometry depends on
 * \tparam ParameterStorageTypes Any number of storage types which describe the parameters
 */
template< typename Kernel, template<class> class Input,
          template<class> class Output, typename... ParameterStorageTypes>
struct DependendGeometry : public UnaryEquation<Kernel, Input<Kernel>, Output<Kernel>>  {
    
    typedef UnaryEquation<Kernel, Input<Kernel>, Output<Kernel>>                   Inherited;
    typedef typename Kernel::Scalar                                                Scalar;
    typedef typename geometry::Geometry<Kernel, ParameterStorageTypes...>::Storage ParameterStorage;

    DependendGeometry() {
        fusion::for_each(m_parameterStorage, detail::Counter<Kernel>(Inherited::m_parameterCount));
    };
    
    //make sure the parameter storage is used, not the geometry one, for initialisation
    virtual void init(LinearSystem<Kernel>& sys) {
#ifdef DCM_DEBUG
        dcm_assert(!Inherited::m_init);
        Inherited::m_init = true;
#endif
        
        //we are a unary equation, we need to make sure that we take our responsibility of 
        //ownership serious
        dcm_assert(Inherited::m_input);
        if(Inherited::hasInputOwnership())
            Inherited::m_input->init(sys);
        
        //mpl trickery to get a sequence counting from 0 to the size of stroage entries
        typedef mpl::range_c<int,0, mpl::size<ParameterStorage>::value> StorageRange;
        
        //now iterate that sequence so we can access all storage elements with knowing the position
        //we are at (that is important to access the correct derivative storage position too)
        mpl::for_each<StorageRange>(detail::Initializer<Kernel, ParameterStorage, Inherited>(sys,
                                    m_parameterStorage, Inherited::m_parameters, 
                                    Inherited::m_derivatives));
        
        //now add the derivatives we take over from the input geometry
        Inherited::m_derivatives.clear();
        for(const auto& param : Inherited::inputEquation()->parameters()) 
            Inherited::m_derivatives.push_back(std::make_pair(typename Inherited::OutputType(), param));
    };
    
    CALCULATE() {
        //ensure the dependend input is calculated correctly
        dcm_assert(Inherited::m_input);        
        if(Inherited::hasInputOwnership())
            Inherited::m_input->execute(); 
                
        //to calculate the real output one needs to know the mathematical equation... this must be done 
        //by the derived class
    };
    
protected:
    ParameterStorage     m_parameterStorage;

    int m_counter = -1; //we need a counter for every calculate, and we do not want a memory allocation 
                        //for every recalculate
};

} //numeric

namespace symbolic {

struct Geometry {
    
    int type;
};

template<typename Kernel, template<class> class G>
struct TypeGeometry : public Geometry {

    typedef G<Kernel> PrimitiveGeometry; 
    
    PrimitiveGeometry& getPrimitveGeometry() {return m_geometry;};
    
    void setPrimitiveGeometry(const G<Kernel>& g) {
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
