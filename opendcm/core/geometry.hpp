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
#include <boost/mpl/transform.hpp>
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

class GeometryIndex {
protected:      
    static int generateIndex(){
        static int idx = -1;
        return ++idx;
    }
};

/**
 * @brief Helper base struct for creating primitive gepometry types
 *
 * Primitive geometry has to fullfill three requirements: the storage based on fusion::vector, the 
 * transform interface and the copy-constructability and assignability. To reduce boilerplate this class
 * is given which handles the fusion::vector creation automatically. To tell the struct
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
 * @note To allow easy access to the stored elements the template function at<Idx> is given, which returnes
 *       the stored data at the given idex.
 *
 * \tparam Kernel The math \ref Kernel in use
 * \tparam StorageTypes Variadic sequence of storage types
 */
template<typename Kernel, typename... StorageTypes>
struct Geometry : GeometryIndex {

    typedef mpl::vector< StorageTypes... >                               StorageSequence;
    typedef typename fusion::result_of::as_vector<StorageSequence>::type Storage;
    
    Geometry& operator=(const Storage& storage) {m_storage = storage; return *this;};
    
    /**
     * @brief Returns the index of the geometry
     */
    static int index() {
        static int idx = generateIndex();
        return idx;
    }
protected:
    Storage m_storage;
    
    template<int Idx>
    auto at()->decltype(fusion::at_c<Idx>(m_storage)) {
        return fusion::at_c<Idx>(m_storage);
    };
};

/**
 * @brief Create a assignable storage for primitive Geometries
 * 
 * Helper function to easily group objects for assignement to an primitive geometry
 */
template<typename... StorageTypes> 
typename fusion::result_of::make_vector<StorageTypes...>::type make_storage(const StorageTypes&... args) {
    return fusion::make_vector(args...);
}

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
    
    using type = Base<Kernel>;
    template<typename K> using  primitive = Base<K>;
};

}//geometry


namespace details {
//helper classes for numeric geometry
template<typename Equation, typename StorageType, typename IdStorageType, bool InitDerivative = true>
struct Initializer {

    typedef typename Equation::KernelType       Kernel;
    typedef typename Kernel::Scalar             Scalar;
    typedef typename Equation::Parameter        Parameter;
    
    numeric::LinearSystem<Kernel>&                  m_system;
    std::vector<Parameter>&                         m_parameters;
    std::vector<typename Equation::DerivativePack>& m_derivatives;
    StorageType&                                    m_storage;
    IdStorageType&                                  m_ids; 
    std::vector<typename Equation::Id>&             m_fixedOutputs;
    
    Initializer(Equation* eqn,
                numeric::LinearSystem<Kernel>& s, 
                StorageType& st, 
                IdStorageType& idstore)                        
        : m_system(s), m_parameters(eqn->m_parameters), m_derivatives(eqn->m_derivatives), m_storage(st), 
          m_fixedOutputs(eqn->m_fixedOutputs), m_ids(idstore) {};

    template<typename T>
    void operator()(T t) {

        //get the parameter
        auto& t1 = fusion::at<T>(m_storage);
        auto& id = fusion::at<T>(m_ids);
        map<T>(t1, id);
    }
    
    template<typename Idx, typename T, typename T2>
    void map(Eigen::MatrixBase<T>& m, Eigen::MatrixBase<T2>& id) {
        
        for(int i=0; i<m.cols(); ++i) {
            for(int j=0; j<m.rows(); ++j) {
                //we need to ensure that only unfixed outputs get a free parameter
                if(!isFixed(id(j,i))) {
                    
                    m_parameters.push_back(Parameter(m_system.mapParameter(), &m(j,i)));
                    m_derivatives.emplace_back(typename Equation::DerivativePack(typename Equation::OutputType(), m_parameters.back()));
                    if(InitDerivative) {
                        auto& t2 = fusion::at<Idx>(m_derivatives.back().first.m_storage);
                        setOne(t2,j,i);
                    }
                }
            }
        }                
    };
    
    template<typename Idx>
    void map(Scalar& s, typename Equation::Id& id) {
        
        if(!isFixed(id)) {
            m_parameters.push_back(Parameter(m_system.mapParameter(), &s));
            m_derivatives.emplace_back(typename Equation::DerivativePack(typename Equation::OutputType(),  m_parameters.back()));
            if(InitDerivative) {
                auto& t2 = fusion::at<Idx>(m_derivatives.back().first.m_storage);
                setOne(t2,0,0);
            }
        }
    };
    
    bool isFixed(typename Equation::Id& id) {
        return std::find(m_fixedOutputs.begin(), m_fixedOutputs.end(), id) != m_fixedOutputs.end();
    };
    
    template<typename T>
    void setOne(Eigen::MatrixBase<T>& m, int row, int col) {
        m(row, col) = 1;
    };
    void setOne(Scalar& s, int row, int col) {
        s = 1;
    };
};

template<typename Kernel>
struct Counter {

    typedef typename Kernel::Scalar Scalar;
    
    unsigned int& count;
    Counter(unsigned int& c) : count(c) {
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

//ID map handling: set everything to undefined
template<typename Id>
struct IdInitalizer {
       
    template<typename T>
    void operator()(Eigen::MatrixBase<T>& val) const {
        val.setConstant(-1);        
    };
    void operator()(Id& s) const {
        s = -1;
    };
};

//ID map handling: find id in given storage
template<typename Id>
struct IdFinder {

    const Id& m_id;
    bool&     m_result;
    
    IdFinder(const Id& id, bool& result) : m_id(id), m_result(result) {};
        
    template<typename T>
    void operator()(Eigen::MatrixBase<T>& val) const {
        if((val.array() == m_id).any())
            m_result = true;
    };
    void operator()(Id& s) const {
        if(s == m_id)
            m_result = true;
    };
};

};//details
    
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
    using Inherited::operator=;

    //mpl trickery to get a sequence counting from 0 to the size of stroage entries
    typedef mpl::range_c<int,0,
            mpl::size<typename Inherited::StorageSequence>::value> StorageRange;
                
    Geometry() {
        Inherited::m_complexity = Complexity::Simple;
        fusion::for_each(Inherited::m_storage, details::Counter<Kernel>(Inherited::m_parameterCount));
    };
    virtual ~Geometry(){};
    
    /** 
     * @brief Returns the maximal parameters this nodes geometry needs
     * Note that this is the static parameter count, meaning its the maximum. It does not consider 
     * any parameter reductions due to unary constraints. 
     */
    static unsigned int staticParameterCount() {
        unsigned int val;
        mpl::for_each<typename Inherited::Storage>(details::Counter<Kernel>(val));
        return val;
    };

    //sometimes it is possible to optimize constraint derivative calculation when we are
    //a normal geometry (where parameret == value). To enable those optimisations we need
    //an indicator if this situation is given
    bool isIndependend() {
        return m_independent;
    };

    //the assumptions made here are only valid for pure geometry, so the derived
    //types must override the initialisaion behaviour
    virtual void init(LinearSystem<Kernel>& sys) override {
#ifdef DCM_DEBUG
        dcm_assert(!Inherited::m_init);
        Inherited::m_init = true;
#endif
        //mpl trickery to get a sequence counting from 0 to the size of stroage entries
        typedef mpl::range_c<int,0,
                mpl::size<typename Inherited::StorageSequence>::value> StorageRange;
        //now iterate that sequence so we can access all storage elements with knowing the position
        //we are at (that is important to access the correct derivative storage position too)
        typedef typename Inherited::OutputIdType::Storage IdStorage;
        mpl::for_each<StorageRange>(details::Initializer<Geometry, typename Inherited::Storage, IdStorage>(
                                        this, sys, Inherited::m_storage, Inherited::m_outputId.m_storage));
        
        //setup the parameter values
        this->storageToParam();
    };
    
    //we actually do not really need to calculate anything, but we need to make sure the mapped 
    //values are moved over to the output
    CALCULATE() {
        this->paramToStorage();
    };
    
protected:
    
    bool m_independent = true;
    
protected:
    int m_counter = -1; //we need a counter for every calculate, and we do not want a memory allocation 
                        //for every recalculate
                           
    using Inherited::m_parameters;
    using Inherited::m_derivatives;
    using Inherited::m_fixedOutputs;
    
    template<typename, typename, typename, bool>
    friend struct details::Initializer;
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
    
private:
    typedef typename mpl::transform<mpl::vector<ParameterStorageTypes...>, details::IdentifierOutputType<mpl::_1>>::type IdVector;
    typedef typename fusion::result_of::as_vector<IdVector>::type ParameterIdStorage;

public:    
    typedef typename Kernel::Scalar         Scalar;
    typedef Equation<Kernel, Base<Kernel>>  Inherited;
    typedef typename geometry::Geometry<Kernel, ParameterStorageTypes...>::Storage ParameterStorage;
    using Inherited::operator=;
    
    ParameterGeometry() {
        Inherited::m_complexity = Complexity::Complex;
        fusion::for_each(m_parameterStorage, details::Counter<Kernel>(Inherited::m_parameterCount));
        fusion::for_each(m_parameterIdStorage, details::IdInitalizer<typename Inherited::Id>());
    };
    virtual ~ParameterGeometry(){};
    
    /** 
     * @brief Returns the maximal parameters this nodes geometry needs
     * Note that this is the static parameter count, meaning its the maximum. It does not consider 
     * any parameter reductions due to unary constraints. 
     */
    static unsigned int staticParameterCount() {
        unsigned int val;
        mpl::for_each<ParameterStorage>(details::Counter<Kernel>(val));
        return val;
    };
    
    //make sure the parameter storage is used, not the geometry one, for initialisation
    virtual void init(LinearSystem<Kernel>& sys) override {
#ifdef DCM_DEBUG
        dcm_assert(!Inherited::m_init);
        Inherited::m_init = true;
#endif
        //mpl trickery to get a sequence counting from 0 to the size of stroage entries
        typedef mpl::range_c<int,0, mpl::size<ParameterStorage>::value> StorageRange;
        
        //now iterate that sequence so we can access all storage elements with knowing the position
        //we are at (that is important to access the correct derivative storage position too)
        mpl::for_each<StorageRange>(details::Initializer<ParameterGeometry, ParameterStorage, ParameterIdStorage, false>(
                                                            this, sys, m_parameterStorage, m_parameterIdStorage));
        
        //setup the parameter values correctly 
        this->storageToParam();
    };
    
    CALCULATE() {
        this->paramToStorage();  
    };
    
    template<int I>
    typename fusion::result_of::at_c<ParameterIdStorage, I>::type parameterIdAt() {
        return fusion::at_c<I>(m_parameterIdStorage);
    };
    
    virtual bool canFixOutput(const typename Inherited::Id& id) override {
        bool result = false;
        fusion::for_each(m_parameterIdStorage, details::IdFinder<typename Inherited::Id>(id, result));
        return result;
    };

protected:
    
    ParameterStorage                                               m_parameterStorage;
    ParameterIdStorage                                             m_parameterIdStorage;
    std::vector<std::pair<typename Inherited::Parameter, Scalar*>> m_paramStorageMap;

    int m_counter = -1; //we need a counter for every calculate, and we do not want a memory allocation 
                        //for every recalculate
                        
    using Inherited::m_parameters;
    using Inherited::m_derivatives;
    using Inherited::m_fixedOutputs;
    
    template<typename, typename, typename, bool>
    friend struct details::Initializer;
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
 * \tparam Output The geometric primitive this numeric geometry represents
 * \tparam ParameterStorageTypes Any number of storage types which describe the parameters
 */
template< typename Kernel, template<class> class Input,
          template<class> class Output, typename... ParameterStorageTypes>
struct DependendGeometry : public UnaryEquation<Kernel, Input<Kernel>, Output<Kernel>>  {
    
private:
    typedef typename mpl::transform<mpl::vector<ParameterStorageTypes...>, details::IdentifierOutputType<mpl::_1>>::type IdVector;
    typedef typename fusion::result_of::as_vector<IdVector>::type ParameterIdStorage;

public:
    typedef UnaryEquation<Kernel, Input<Kernel>, Output<Kernel>>                   Inherited;
    typedef typename Kernel::Scalar                                                Scalar;
    typedef typename geometry::Geometry<Kernel, ParameterStorageTypes...>::Storage ParameterStorage;
    using Inherited::operator=;

    DependendGeometry() {
        Inherited::m_complexity = Complexity::Complex;
        fusion::for_each(m_parameterStorage, details::Counter<Kernel>(Inherited::m_parameterCount));
        fusion::for_each(m_parameterIdStorage, details::IdInitalizer<typename Inherited::Id>());       
    };
    virtual ~DependendGeometry(){};
    
    /** 
     * @brief Returns the maximal parameters this nodes geometry needs
     * Note that this is the static parameter count, meaning its the maximum. It does not consider 
     * any parameter reductions due to unary constraints. 
     */
    static unsigned int staticParameterCount() {
        unsigned int val;
        mpl::for_each<ParameterStorage>(details::Counter<Kernel>(val));
        return val;
    };
    
    //make sure the parameter storage is used, not the geometry one, for initialisation
    virtual void init(LinearSystem<Kernel>& sys) override {
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
        mpl::for_each<StorageRange>(details::Initializer<DependendGeometry, ParameterStorage, ParameterIdStorage, false>(
                                        this, sys, m_parameterStorage, m_parameterIdStorage));
        
        //now add the derivatives we take over from the input geometry
        Inherited::m_derivatives.clear();
        for(const auto& param : Inherited::inputEquation()->parameters()) 
            Inherited::m_derivatives.push_back(std::make_pair(typename Inherited::OutputType(), param));
        
        //make sure the parameters have the correct values
        this->storageToParam();
    };
    
    CALCULATE() {
        //ensure the dependend input is calculated correctly
        dcm_assert(Inherited::m_input);        
        if(Inherited::hasInputOwnership())
            Inherited::m_input->execute(); 
                
        //to calculate the real output one needs to know the mathematical equation... this must be done 
        //by the derived class
    };
    
    template<int I>
    typename fusion::result_of::at_c<ParameterIdStorage, I>::type parameterIdAt() {
        return fusion::at_c<I>(m_parameterIdStorage);
    };
    
    virtual bool canFixOutput(const typename Inherited::Id& id) override {
        bool result = false;
        fusion::for_each(m_parameterIdStorage, details::IdFinder<typename Inherited::Id>(id, result));
        return result;
    };
    
protected:
    
    ParameterStorage                                               m_parameterStorage;
    ParameterIdStorage                                             m_parameterIdStorage;
    std::vector<std::pair<typename Inherited::Parameter, Scalar*>> m_paramStorageMap;

    int m_counter = -1; //we need a counter for every calculate, and we do not want a memory allocation 
                        //for every recalculate
                        
    using Inherited::m_parameters;
    using Inherited::m_derivatives;
    using Inherited::m_fixedOutputs;
    
    template<typename, typename, typename, bool>
    friend struct details::Initializer;
};

} //numeric

namespace symbolic {

struct Geometry {
    
    int  getType() {return m_type;};
    
protected:
    int m_type;
};

template<typename Primitive>
struct TypeGeometry : public Geometry {

    void       setPrimitive(const Primitive& g) {
        m_geometry = g;
        m_type = g.index();
    }
    Primitive& getPrimitve() {return m_geometry;};
    
protected:   
    Primitive m_geometry;
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
