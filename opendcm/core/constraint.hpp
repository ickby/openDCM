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
    GNU Lesser General Public License for more detemplate tails.

    You should have received a copy of the GNU Lesser General Public License along
    with this library; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/

#ifndef DCM_CONSTRAINT_H
#define DCM_CONSTRAINT_H

#include <boost/utility/enable_if.hpp>
#include <boost/mpl/contains.hpp>
#include <boost/fusion/include/vector.hpp>
#include <boost/fusion/include/find.hpp>
#include <boost/fusion/include/at.hpp>
#include <boost/exception/errinfo_errno.hpp>

#include <iostream>
#include <type_traits>
#include "geometry.hpp"
#include "defines.hpp"

namespace mpl = boost::mpl;
namespace fusion = boost::fusion;

namespace dcm { 
    
//the possible directions
enum class Orientations { Parallel, Equal, Opposite, Perpendicular };

//the possible solution spaces
enum class SolutionSpaces {Bidirectional, Positiv_directional, Negative_directional};

//the fixable entities
enum class Fixables {pointX, pointY, pointZ, directionX, directionY, directionZ, radius};

/**
 * @brief Constraint primitives handling
 *
 * Equivalent to the geometric primitives there is a definition for constraint primitives. A constraint 
 * primitive is a constraint type like distance, primitive in this context means that it holds only its type
 * and all needed information to fully define it, called option. For a distance constraint this is for 
 * example the nuemeric distance value. A primitive is used to store and access the type and its data in a
 * compfortable and user friendly way.
 * There are a few requirements for primitive constraints. First they must be assignalble and copyconstructible
 * while preserving their option values. Second all options need to be stored in a  fusion vector called
 * m_storage. Third it must be possible to provide options via operator() for one or many options at once 
 * and furthermore through assigning with operator=, also for single options or multiple ones via initializer
 * lists. Fourth the function geometryCount() must return the amount of geometries needed by this constrain. 
 * The last requirement regards the default values of the options. It is important to have the possibility
 * to reset all options to default during the objects lifetime. Therefore the function setDefault() must exist.
 * An example primitive constraint which holds all requirements could look like this:
 * @code
 * struct TestConstraint {
 *   
 *  fucion::vector<int, char> m_sequence; //requirement 2
 *  
 *  //requirement 3
 *  void operator()(int);
 *  void operator()(char);
 *  void operator()(int, char);
 *  void operator=(int);
 *  void operator=(char);
 *  void operator=(int, char);
 * 
 *  int geometryCount(); //requirement 4;
 *  void setDefault(); //requirement 5;
 * };
 * @endcode
 * 
 * Requirement 1 is fullfilled by the standart compiler generated assignement operator and copy constructor. This
 * also holds for the assignemnd of multiple options via initializer lists.
 * 
 * \Note The multiassignement via operator() must not support arbitrary argument ordering. It is sufficient to 
 * allow the parameter only in the same order as the options are specified.
 */
namespace constraint {  

template<int Arity>
class ConstraintIndex {
protected:      
    static int generateIndex(){
        static int idx = -1;
        return ++idx;
    }
};

/**
 * @brief Herlper class to create primitive constraints
 * 
 * The primitive constraint type concept formulates a few requirements on the used types. To achieve a full
 * compatibility with the concept a certain boilerplate is needed. To circumvent this boilerplate this class
 * is given. It accepts the option types and the geometry count for the constraint as template arguments and
 * then provides all needed storages and assignment operators as well as all needed functions. 
 * This only works for option types which are default constructible. The only second restriction is that the
 * option types need to distuinguishable. To use a option type twice is not allowd, for example two times int.
 * Otherwise the option assignement of single options would be ambigious. 
 * The default value requirement is achieved through the constructor which acceppts the default values. If 
 * the standart constructor is used then the default constructed values of the option types are used as 
 * default values. Note that you must provide the constructors also in your derived class.
 * An example constraint can look like this:
 * @code
 * struct TetstConstraint : public constraint::Constraint<TestConstraint, 2, int, char> {
 *    using Constraint::operator=;
 *    TestConstraint() {};
 *    TestConstraint(const int& i, const char& c) : Constraint(i,c) {};
 * };
 * TestConstraint test(1,'a');
 * @endcode
 * @note To allow easy access to the stored options the template function at<Idx> is given, which returnes
 *       the stored data at the given idex. 
 * \tparam Derived The type of the derived constraint 
 * \tparam GCount The integer defining the amount of geometies needed for this constraint
 * \tparam OptionTypes a variadic sequence of copy constructable types which describe the stored options
*/    
template<typename Derived, int GCount, typename ...OptionTypes>
struct Constraint : public ConstraintIndex<GCount> {

    const static int Arity = GCount;
    typedef typename fusion::vector<OptionTypes...> Options;
    typedef typename mpl::at_c<Options, 0>::type    PimaryOptionType;

    Constraint() {};
    
    Constraint(const OptionTypes&... defaults) : m_defaults(defaults...) {
        m_storage = m_defaults;
    };
    
    Derived& derived() {return *static_cast<Derived*>(this);};
    const Derived& derived() const {return *static_cast<const Derived*>(this);};
    
    //Copy assign option. We provide this to allow the automatic copy constructor to be generated. 
    //This gives 2 operators for the price of implementing one, for this base and all derived classes.
    //As we have only one parameter in this function it would be ambigious with the single option 
    //operator(), hence we need to activate/deactivate both according to the provided type.
    template<typename T>
    typename boost::disable_if<mpl::contains<Options, T>, Derived&>::type operator()(const T& c) {
        m_storage = c.m_storage;
        return derived();
    };

    //Set a single option. As we have only one parameter in this function it would be ambigious with the 
    //copy operator(), hence we need to activate/deactivate both according to the provided type.
    template<typename T>
    typename boost::enable_if<mpl::contains<Options, T>, Derived&>::type operator()(const T& val) {
        BOOST_MPL_ASSERT((mpl::contains<Options, T>));
        *fusion::find<T>(m_storage) = val;
        return derived();
    };
    
    //set multiple options at once.
    Derived& operator()(const OptionTypes&... val) {
        m_storage = Options(val...);
        return derived();
    };

    //Assign option. Disable it to avoid confusion with the copy assignement operator
    template<typename T>
    typename boost::enable_if<mpl::contains<Options, T>, Derived&>::type operator=(const T& val) {
        return operator()(val);
    };
    
    //explicity give copy assignement as it otherwise would be treated as deleted due to the assignend operator
    //given above. When this one is provided the derived class wil also generate one.
    Constraint& operator=(const Constraint& val) {
        m_storage = val.m_storage;
        return *this;
    };
   
    //set default option values, neeeded for repedability and to prevent unexpected behaviour
    void setDefault() {
        m_storage = m_defaults;
    };
    
    int geometryCount() {
        return GCount;
    };
    
    //access a option by index 
    template<int idx>
    typename fusion::result_of::at_c<Options, idx>::type getOption() {
        return fusion::at_c<idx>(m_storage);
    };
    
    const Options& getOptions() {
        return m_storage;
    };
    
    /**
     * @brief Returns the index of the geometry
     */
    static int index() {
        static int idx = ConstraintIndex<Arity>::generateIndex();
        return idx;
    }
       
protected:
    Options       m_storage;
    const Options m_defaults;
    
    template<int Idx>
    auto at()->decltype(fusion::at_c<Idx>(m_storage)) {
        return fusion::at_c<Idx>(m_storage);
    };
};

}//constraint
    
    
    
namespace symbolic {
    
struct Constraint {
    
    void setType(int id) { type = id;};
    int  getType() {return type;};
    
    void setArity(int a) {arity = a;};
    int  getArity() {return arity;};
    
protected:
    int type, arity;
};

template<typename Primitive>
struct TypeConstraint : public Constraint {

    void       setPrimitive(const Primitive& c) {m_constraint = c;};    
    Primitive& getPrimitive() {return m_constraint;};    
    
protected:   
    Primitive m_constraint;
};

struct ConstraintProperty {
    typedef Constraint* type;
    struct default_value {
        Constraint* operator()() {
            return nullptr;
        };
    };
    struct change_tracking{};
};

struct ConstraintListProperty {
    typedef std::vector<Constraint*> type;
    struct change_tracking{};
};
    
}//symbolic    
    
namespace numeric {
 
/**
 * @brief Base class to unify derivation of parent classes: Single geometry constraints
 * 
 * To ease the specification of the parent classes for specilized numeric constraint classes as well 
 * as simplifying the needed types this class is given. It encapsulates all needed inheritance and 
 * typedefs. Note that the equation exposes a single scalar as result. This is due to the fact, that 
 * a constraint equation is always a error function for the numeric solver. Hence this class also 
 * ensures that the Calculatable newResidualCount() returns one. 
 */
template<typename Kernel, typename PC, typename PG>
struct UnaryConstraintBase : public UnaryEquation<Kernel, PG, typename Kernel::Scalar>, 
                        public PC  {
    
        typedef UnaryEquation<Kernel, PG, typename Kernel::Scalar> Equation;
        
        typedef typename Kernel::Scalar                  Scalar;
        typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Vector;
        typedef PG                                       Geometry;
        typedef Geometry                                 Derivative;
        typedef typename Equation::InputEqn              InputEquation;
        
        using PC::operator=;
        
        void assign(const PC& pc) {
            operator=(pc);
        };
        
        UnaryConstraintBase() {
            Equation::m_residualCount = 1;
        };
};

/**
 * @brief Base class to unify derivation of parent classes: Double geometry equations
 * 
 * To ease the specification of the parent classes for specilized numeric constraint classes as well 
 * as simplifying the needed types this class is given. It encapsulates all needed inheritance and 
 * typedefs. Note that the equation exposes a single scalar as result. This is due to the fact, that 
 * a constraint equation is always a error function for the numeric solver. Hence this class also 
 * ensures that the Calculatable newResidualCount() returns one. 
 */
template<typename Kernel, typename PC, typename PG1, typename PG2>
struct BinaryConstraintBase : public BinaryEquation<Kernel, PG1, PG2, typename Kernel::Scalar>, 
                        public PC  {
    
        typedef BinaryEquation<Kernel, PG1, PG2, typename Kernel::Scalar> Equation;
        
        typedef typename Kernel::Scalar                  Scalar;
        typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Vector;
        typedef PG1                                      Geometry1;
        typedef Geometry1                                Derivative1;
        typedef PG2                                      Geometry2;
        typedef Geometry2                                Derivative2;
        
        using PC::operator=;
        
        void assign(const PC& pc) {
            operator=(pc);
        };
        
        BinaryConstraintBase() {
            Equation::m_residualCount = 1;
        };
};
    
/**
 * @brief Class for numeric evaluation of primitive constraints
 * 
 * Primitive constraints hold only their constraint type and the options to fully define their behaviour. 
 * As numeric evaluation of constraints may be needed in the solving process this information is not enough.
 * The governing equations for the primitive constraints and their various possible geometry combinations 
 * are needed. This and all derived classes provide a conviniet interface to handle all required equations.
 * 
 * This class provides a way of specifying equations for numeric evaluation of constraints. As
 * every primitive constraint can be defined for multiple geometry combinations it must allow to specify 
 * the equations for those different contexts. To achieve this a numeric constraint is a template class,
 * with the primitive constraint and the geometries as template parameters. This allows to specialize the
 * class for every possible combination.
 *
 * To fully evaluate a equation in the numeric solving process multiple functions are needed. For one the 
 * error function, but also the functions to evaluate the error functions derivative. The derivative can
 * in its simplest form be calculated for every single parameter in the geometries. However, even if 
 * always correct, this is often more work than needed. Often, a parameter is directly used in the given
 * equation and hence only one equation parameter will depend on one geometry parameter. This allows to
 * dramatically reduce the amount of calculations needed for evaluating the derivative. Therefore this 
 * class offers a way to specify both, the full evaluation and the optimized one.
 *
 * In order to give the equations for a certain combination one provides the following code:
 * 
 * \code{.cpp}
 * 
 * namespace dcm { namespace numeric {
 * 
 * template<typename Kernel>
 * struct BinaryConstraint<Kernel, dcm::Distance, dcm::Point3, dcm::Point3> : public BinaryConstraintBase<Kernel, dcm::Distance, TPoint3, TPoint3>  {
 *
 *      Scalar calculateError(const Geometry1& g1, const Geometry2& g2) {};
 *      Scalar calculateGradientFirst(const Geometry1& g1, const Geometry2& g2, const Derivative1& dg1) {};
 *      Scalar calculateGradientSecond(const Geometry1& g1, const Geometry2& g2, const Derivative2& dg2) {};
 *      Vector calculateGradientFirstComplete(const Geometry1& g1, const Geometry2& g2) {};
 *      Vector calculateGradientSecondComplete(const Geometry1& g1, const Geometry2& g2) {};
 * }
 * }}
 * \endcode
 *  
 * \note The numeric constraint is derived from the primitive constraint, therefore one can use
 * the primitive constraint functions to access the constraint options.
 * 
 * \tparam Kernel the math kernel in use
 * \tparam PC  the primitive constraint in use
 * \tparam PG1 the first primitive geometry the equation is defined for
 * \tparam PG2 the second primitive geometry the equation is defined for
 */
template<typename Kernel, typename PC, typename PG1, typename PG2>
struct BinaryConstraint : public BinaryConstraintBase<Kernel, PC, PG1, PG2> {
           
        typedef BinaryConstraintBase<Kernel, PC, PG1, PG2> Inherited;
        typedef typename Kernel::Scalar                    Scalar;
        typedef typename Inherited::Vector                 Vector;
        typedef typename Inherited::Geometry1              Geometry1;
        typedef typename Inherited::Derivative1            Derivative1;
        typedef typename Inherited::Geometry2              Geometry2;
        typedef typename Inherited::Derivative2            Derivative2;
        
        typedef PC  InheritedConstraint;
        typedef BinaryEquation<Kernel, Geometry1, Geometry2, Scalar> InheritedEquation;
        
        BinaryConstraint() {
            throw creation_error() <<  boost::errinfo_errno(24) << error_message("Constraint is not supported for given geometry types");
        };
        
        Scalar calculateError(const Geometry1& g1, const Geometry2& g2) {
            throw creation_error() <<  boost::errinfo_errno(24) << error_message("Constraint is not supported for given geometry types");
        };

        Scalar calculateGradientFirst(const Geometry1& g1, const Geometry2& g2, const Derivative1& dg1) {
            throw creation_error() <<  boost::errinfo_errno(24) << error_message("Constraint is not supported for given geometry types");
        };

        Scalar calculateGradientSecond(const Geometry1& g1, const Geometry2& g2, const Derivative2& dg2) {
            throw creation_error() <<  boost::errinfo_errno(24) << error_message("Constraint is not supported for given geometry types");
        };

        Vector calculateGradientFirstComplete(const Geometry1& g1, const Geometry2& g2) {
            throw creation_error() <<  boost::errinfo_errno(24) << error_message("Constraint is not supported for given geometry types");
        };

        Vector calculateGradientSecondComplete(const Geometry1& g1, const Geometry2& g2) {
            throw creation_error() <<  boost::errinfo_errno(24) << error_message("Constraint is not supported for given geometry types");
        };
};

template<typename Kernel, typename PC, typename PG>
struct UnaryConstraint : public UnaryConstraintBase<Kernel, PC, PG> {
           
        typedef UnaryConstraintBase<Kernel, PC, PG>      Inherited;
        typedef typename Kernel::Scalar                  Scalar;
        typedef typename Inherited::Vector               Vector;
        typedef typename Inherited::Geometry             Geometry;
        typedef typename Inherited::Derivative           Derivative;

        typedef PC                                       Constraint;
        typedef typename Inherited::InputEquation        InputEquation;
        
        UnaryConstraint() {
            throw creation_error() <<  boost::errinfo_errno(24) << error_message("Constraint is not supported for given geometry types");
        };
        
        bool applyToEquation(std::shared_ptr<InputEquation> eqn) {
            throw creation_error() <<  boost::errinfo_errno(24) << error_message("Constraint is not supported for given geometry types");
        }
        
        Scalar calculateError(const Geometry& g1) {
            throw creation_error() <<  boost::errinfo_errno(24) << error_message("Constraint is not supported for given geometry types");
        };

        Scalar calculateGradient(const Geometry& g, const Derivative& dg) {
            throw creation_error() <<  boost::errinfo_errno(24) << error_message("Constraint is not supported for given geometry types");
        };

        Vector calculateGradientComplete(const Geometry& g) {
            throw creation_error() <<  boost::errinfo_errno(24) << error_message("Constraint is not supported for given geometry types");
        };
};
    
/**
* @brief Numeric handling of error functions
* 
* As all error functions share a certain kind of structure, this class is used to provide a single 
* implementation for reused functionality. This involves the initialisation of the equation with 
* allocating the residual and the derivatives in the solver. This class provides the storage and 
* access points for the result and derivatives. It further provides higher level functions to calculate 
* individual parts of the equation.
* 
*/    
template<typename Kernel, typename PC, typename PG1, typename PG2>
struct BinaryConstraintEquationBase : public numeric::BinaryConstraint<Kernel, PC, PG1, PG2> {
   
    typedef BinaryConstraint<Kernel, PC, PG1, PG2>      Inherited;
    typedef VectorEntry<Kernel>                         Residual;
    typedef MatrixEntry<Kernel>                         Derivative;
    
    //type to hold geometric derivative together with the correct position for the jacobi entry
    typedef std::pair<typename numeric::Equation<Kernel, PG1>::Derivative*, Derivative> Derivative1Pack;
    typedef std::pair<typename numeric::Equation<Kernel, PG2>::Derivative*, Derivative> Derivative2Pack;
     
    virtual void init(LinearSystem<Kernel>& sys) {
#ifdef DCM_DEBUG
        dcm_assert(!m_init);
        dcm_assert(Inherited::firstInputEquation() && Inherited::firstInputEquation()->isInitialized());
        dcm_assert(Inherited::secondInputEquation() && Inherited::secondInputEquation()->isInitialized());
        m_init = true;
#endif
        //setup the residual first to see in which row we are working with this constraint
        residual = sys.mapResidual();
            
        //Setup the correct jacobi entry for the individual parameter
        for(auto& der : Inherited::firstInputEquation()->derivatives())  
            g1_derivatives.push_back({&der.first, sys.mapJacobi(residual.Index, der.second.getEntry().Index)});
    
        for(auto& der : Inherited::secondInputEquation()->derivatives())
            g2_derivatives.push_back({&der.first, sys.mapJacobi(residual.Index, der.second.getEntry().Index)});
    };
    
#ifdef DCM_TESTING
    typename Kernel::Scalar getResidual() {
        return *residual.Value;
    };
    
    std::vector<Derivative1Pack> getFirstDerivatives() {
        return g1_derivatives;
    }
    
    std::vector<Derivative2Pack> getSecondDerivatives() {
        return g2_derivatives;
    }
#endif
    
protected:
    
    void firstAsSimplified() {
         
        auto result1 = Inherited::calculateGradientFirstComplete(Inherited::firstInput(), 
                                                                 Inherited::secondInput());
        dcm_assert(result1.rows() == g1_derivatives.size());
        int i = 0;
        for(Derivative1Pack& der : g1_derivatives) {
            *(der.second.Value) = result1(i);
            ++i;
        }
    };
    
    void firstAsComplex() {
      
        for(Derivative1Pack& der : g1_derivatives) 
            *(der.second.Value) = Inherited::calculateGradientFirst(Inherited::firstInput(),
                                                                    Inherited::secondInput(), *der.first);
    };
    
    void secondAsSimplified() {
        
        auto result2 = Inherited::calculateGradientSecondComplete(Inherited::firstInput(), 
                                                                  Inherited::secondInput());
        dcm_assert(result2.rows() == g2_derivatives.size());
        int i = 0;
        for(Derivative2Pack& der : g2_derivatives) {
            *(der.second.Value) = result2(i);
            ++i;
        }
    };
    
    void secondAsComplex() {
        
        for(Derivative2Pack& der : g2_derivatives) 
            *(der.second.Value) = Inherited::calculateGradientSecond(Inherited::firstInput(),
                                                                     Inherited::secondInput(), *der.first);
    };

#ifdef DCM_DEBUG
    bool m_init = false;
#endif
    Residual                        residual;
    std::vector<Derivative1Pack>    g1_derivatives;
    std::vector<Derivative2Pack>    g2_derivatives;
    
};

    
/**
* @brief Numeric handling of error functions for single geometry constraints
* 
* As all error functions share a certain kind of structure, this class is used to provide a single 
* implementation for reused functionality. This involves the initialisation of the equation with 
* allocating the residual and the derivatives in the solver. This class provides the storage and 
* access points for the result and derivatives. It further provides higher level functions to calculate 
* individual parts of the equation.
* 
*/    
template<typename Kernel, typename PC, typename PG>
struct UnaryConstraintEquationBase : public numeric::UnaryConstraint<Kernel, PC, PG> {
   
    typedef UnaryConstraint<Kernel, PC, PG>  Inherited;
    typedef VectorEntry<Kernel>              Residual;
    typedef MatrixEntry<Kernel>              Derivative;
    
    //type to hold geometric derivative together with the correct position for the jacobi entry
    typedef std::pair<typename numeric::Equation<Kernel, PG>::Derivative*, Derivative> DerivativePack;
     
    virtual void init(LinearSystem<Kernel>& sys) {
#ifdef DCM_DEBUG
        dcm_assert(!m_init);
        dcm_assert(Inherited::inputEquation() && Inherited::inputEquation()->isInitialized());
        m_init = true;
#endif
        //setup the residual first to see in which row we are working with this constraint
        residual = sys.mapResidual();
            
        //Setup the correct jacobi entry for the individual parameter
        for(auto& der : Inherited::inputEquation()->derivatives())  
            derivatives.push_back({&der.first, sys.mapJacobi(residual.Index, der.second.getEntry().Index)});
    };
    
#ifdef DCM_TESTING
    typename Kernel::Scalar getResidual() {
        return *residual.Value;
    };
    
    std::vector<DerivativePack> getDerivatives() {
        return derivatives;
    }
#endif
    
protected:
    
    void asSimplified() {
         
        auto result1 = Inherited::calculateGradientComplete(Inherited::input());
        dcm_assert(result1.rows() == derivatives.size());
        int i = 0;
        for(DerivativePack& der : derivatives) {
            *(der.second.Value) = result1(i);
            ++i;
        }
    };
    
    void asComplex() {
      
        for(DerivativePack& der : derivatives) 
            *(der.second.Value) = Inherited::calculateGradient(Inherited::input(), *der.first);
    };

#ifdef DCM_DEBUG
    bool m_init = false;
#endif
    Residual                        residual;
    std::vector<DerivativePack>     derivatives;   
};

/**
 * @brief Error function evaluation for two simple inputs
 * 
 * If both input equations are simple, meaning they are not an InputEquation themself, one can use 
 * the simplified functions for both inputs.  This class provides the needed claculation functionality
 * for the equation.
 */
template<typename Kernel, typename PC, typename PG1, typename PG2>
struct BinaryConstraintSimplifiedEquation : public BinaryConstraintEquationBase<Kernel, PC, PG1, PG2> {
    
    typedef BinaryConstraintEquationBase<Kernel, PC, PG1, PG2> Inherited;
    
    CALCULATE() {
        *Inherited::residual.Value = Inherited::calculateError(Inherited::firstInput(), Inherited::secondInput());
        Inherited::firstAsSimplified();
        Inherited::secondAsSimplified();
    };
};

/**
 * @brief Error function evaluation for two complex inputs
 * 
 * If both input equations are complex, meaning they are  InputEquation themself, one must use 
 * the complex functions for both inputs.  This class provides the needed claculation functionality
 * for the equation.
 */
template<typename Kernel, typename PC, typename PG1, typename PG2>
struct BinaryConstraintComplexEquation : BinaryConstraintEquationBase<Kernel, PC, PG1, PG2> {
  
    typedef BinaryConstraintEquationBase<Kernel, PC, PG1, PG2> Inherited;
    
    CALCULATE() {
        *Inherited::residual.Value = Inherited::calculateError(Inherited::firstInput(), Inherited::secondInput());
        Inherited::firstAsComplex();
        Inherited::secondAsComplex();
    };
};

/**
 * @brief Error function evaluation for two complex inputs
 * 
 * If the first input equations is simple, the second complex, one must use 
 * the simple and complex functions respectivly.  This class provides the needed claculation functionality
 * for the equation.
 */
template<typename Kernel, typename PC, typename PG1, typename PG2>
struct BinaryConstraintSimplifiedComplexEquation : BinaryConstraintEquationBase<Kernel, PC, PG1, PG2> {
    
    typedef BinaryConstraintEquationBase<Kernel, PC, PG1, PG2> Inherited;
    
    CALCULATE() {
        *Inherited::residual.Value = Inherited::calculateError(Inherited::firstInput(), Inherited::secondInput());
        Inherited::firstAsSimplified();
        Inherited::secondAsComplex();
    };
};

/**
 * @brief Error function evaluation for two complex inputs
 * 
 * If the second input equations is simple, the first complex, one must use 
 * the simple and complex functions respectivly.  This class provides the needed claculation functionality
 * for the equation.
 */
template<typename Kernel, typename PC, typename PG1, typename PG2>
struct BinaryConstraintComplexSimplifiedEquation : BinaryConstraintEquationBase<Kernel, PC, PG1, PG2> {
  
    typedef BinaryConstraintEquationBase<Kernel, PC, PG1, PG2> Inherited;
   
    CALCULATE() {
        *Inherited::residual.Value = Inherited::calculateError(Inherited::firstInput(), Inherited::secondInput());
        Inherited::firstAsComplex();
        Inherited::secondAsSimplified();
    };
};

/**
 * @brief Error function evaluation for two simple inputs
 * 
 * If both input equations are simple, meaning they are not an InputEquation themself, one can use 
 * the simplified functions for both inputs.  This class provides the needed claculation functionality
 * for the equation.
 */
template<typename Kernel, typename PC, typename PG>
struct UnaryConstraintSimplifiedEquation : public UnaryConstraintEquationBase<Kernel, PC, PG> {
    
    typedef UnaryConstraintEquationBase<Kernel, PC, PG> Inherited;
    
    CALCULATE() {
        *Inherited::residual.Value = Inherited::calculateError(Inherited::input());
        Inherited::asSimplified();
    };
};

/**
 * @brief Error function evaluation for two complex inputs
 * 
 * If both input equations are complex, meaning they are  InputEquation themself, one must use 
 * the complex functions for both inputs.  This class provides the needed claculation functionality
 * for the equation.
 */
template<typename Kernel, typename PC, typename PG>
struct UnaryConstraintComplexEquation : UnaryConstraintEquationBase<Kernel, PC, PG> {
  
    typedef UnaryConstraintEquationBase<Kernel, PC, PG> Inherited;
    
    CALCULATE() {
        *Inherited::residual.Value = Inherited::calculateError(Inherited::input());
        Inherited::asComplex();
    };
};


template<typename Kernel>
struct BinaryConstraintEquationGenerator {
    
     typedef shedule::FlowGraph::Node  FlowNode;
     typedef std::shared_ptr<numeric::Calculatable<Kernel>>  Equation;
    
     virtual ~BinaryConstraintEquationGenerator()  {};
     
     virtual Equation buildEquation(Equation g1, 
                                    Equation g2, 
                                    symbolic::Constraint* c) const = 0;
                                              
    virtual std::pair<Equation, FlowNode>
    buildEquationNode(Equation g1, 
                      Equation g2, 
                      symbolic::Constraint* c,
                      shedule::FlowGraph& flowgraph) const = 0;
                                              
};

template<typename Kernel>
struct UnaryConstraintEquationGenerator {
    
     typedef shedule::FlowGraph::Node  FlowNode;
     typedef std::shared_ptr<numeric::Calculatable<Kernel>>  Equation;
    
     virtual ~UnaryConstraintEquationGenerator()  {};
    
     virtual bool applyToEquation(Equation g, symbolic::Constraint* c) const = 0;
     virtual Equation buildEquation(Equation g, symbolic::Constraint* c) const = 0;
                                              
     virtual std::pair<Equation, FlowNode>
     buildEquationNode(Equation g,
                      symbolic::Constraint* c,
                      shedule::FlowGraph& flowgraph) const = 0;
                                              
};

template<typename Kernel, typename PC, typename PG1, typename PG2>
struct TypedBinaryConstraintEquationGenerator : public BinaryConstraintEquationGenerator<Kernel> {

    typedef typename BinaryConstraintEquationGenerator<Kernel>::Equation Equation;
    typedef typename BinaryConstraintEquationGenerator<Kernel>::FlowNode FlowNode;
    
    virtual ~TypedBinaryConstraintEquationGenerator()  {};
    
    virtual Equation buildEquation(Equation g1, 
                                   Equation g2, 
                                   symbolic::Constraint* c) const override {
        
        auto tg1 = std::static_pointer_cast<numeric::Equation<Kernel, PG1>>(g1);
        auto tg2 = std::static_pointer_cast<numeric::Equation<Kernel, PG2>>(g2);
        auto& pc  = static_cast<symbolic::TypeConstraint<PC>*>(c)->getPrimitive();
        
        if(tg1->getComplexity() != Complexity::Complex && tg2->getComplexity() != Complexity::Complex) {
            auto equation = std::make_shared<BinaryConstraintSimplifiedEquation<Kernel, PC, PG1, PG2>>(); 
            equation->setInputEquations(tg1, tg2);
            equation->assign(pc);            
            return equation;
        }
        else if(tg1->getComplexity() != Complexity::Complex && tg2->getComplexity() == Complexity::Complex) {
            auto equation = std::make_shared<BinaryConstraintSimplifiedComplexEquation<Kernel, PC, PG1, PG2>>();
            equation->setInputEquations(tg1, tg2);            
            equation->assign(pc);
            return equation;
        }
        else if(tg1->getComplexity() == Complexity::Complex && tg2->getComplexity() != Complexity::Complex) {
            auto equation = std::make_shared<BinaryConstraintComplexSimplifiedEquation<Kernel, PC, PG1, PG2>>();
            equation->setInputEquations(tg1, tg2);            
            equation->assign(pc);
            return equation;
        }
        else if(tg1->getComplexity() == Complexity::Complex && tg2->getComplexity() == Complexity::Complex) {
            auto equation = std::make_shared<BinaryConstraintComplexEquation<Kernel, PC, PG1, PG2>>();
            equation->setInputEquations(tg1, tg2);            
            equation->assign(pc);
            return equation;
        }
        dcm_assert(false);
        return Equation();
    };
    
    virtual std::pair<Equation, FlowNode>
    buildEquationNode(Equation g1, Equation g2, symbolic::Constraint* c,
                      shedule::FlowGraph& flowgraph) const override {
       
        auto  tg1 = std::static_pointer_cast<numeric::Equation<Kernel, PG1>>(g1);
        auto  tg2 = std::static_pointer_cast<numeric::Equation<Kernel, PG2>>(g2);
        auto& pc  = static_cast<symbolic::TypeConstraint<PC>*>(c)->getPrimitive();
        
        if(tg1->getComplexity() != Complexity::Complex && tg2->getComplexity() != Complexity::Complex) {
            auto equation = std::make_shared<BinaryConstraintSimplifiedEquation<Kernel, PC, PG1, PG2>>(); 
            equation->setInputEquations(tg1, tg2);            
            equation->assign(pc);
            return std::make_pair(equation, flowgraph.newActionNode([=](const shedule::FlowGraph::ContinueMessage& m){
                equation->calculate();
            }));
        }
        else if(tg1->getComplexity() != Complexity::Complex && tg2->getComplexity() == Complexity::Complex) {
            auto equation = std::make_shared<BinaryConstraintSimplifiedComplexEquation<Kernel, PC, PG1, PG2>>();
            equation->setInputEquations(tg1, tg2);     
            equation->assign(pc);
            return std::make_pair(equation, flowgraph.newActionNode([=](const shedule::FlowGraph::ContinueMessage& m){
                equation->calculate();
            }));
        }
        else if(tg1->getComplexity() == Complexity::Complex && tg2->getComplexity() != Complexity::Complex) {
            auto equation = std::make_shared<BinaryConstraintComplexSimplifiedEquation<Kernel, PC, PG1, PG2>>();
            equation->setInputEquations(tg1, tg2);  
            equation->assign(pc);
            return std::make_pair(equation, flowgraph.newActionNode([=](const shedule::FlowGraph::ContinueMessage& m){
                equation->calculate();
            }));
        }
        else {
            auto equation = std::make_shared<BinaryConstraintComplexEquation<Kernel, PC, PG1, PG2>>();
            equation->setInputEquations(tg1, tg2);  
            equation->assign(pc);
            return std::make_pair(equation, flowgraph.newActionNode([=](const shedule::FlowGraph::ContinueMessage& m){
                equation->calculate();
            }));
        }
    };
};


template<typename Kernel, typename PC, typename PG>
struct TypedUnaryConstraintEquationGenerator : public UnaryConstraintEquationGenerator<Kernel> {

    typedef typename UnaryConstraintEquationGenerator<Kernel>::Equation Equation;
    typedef typename UnaryConstraintEquationGenerator<Kernel>::FlowNode FlowNode;
    
    virtual ~TypedUnaryConstraintEquationGenerator()  {};
    
    virtual bool applyToEquation(Equation g, symbolic::Constraint* c) const  override  {
    
        auto tg = std::static_pointer_cast<numeric::Equation<Kernel, PG>>(g);
        auto& pc  = static_cast<symbolic::TypeConstraint<PC>*>(c)->getPrimitive();
        UnaryConstraint<Kernel, PC, PG> cons;
        cons.assign(pc);
        return cons.applyToEquation(tg);
    };
    
    virtual Equation buildEquation(Equation g, symbolic::Constraint* c) const override {
        
        auto tg = std::static_pointer_cast<numeric::Equation<Kernel, PG>>(g);
        auto& pc  = static_cast<symbolic::TypeConstraint<PC>*>(c)->getPrimitive();
        
        if(tg->getComplexity() != Complexity::Complex) {
            auto equation = std::make_shared<UnaryConstraintSimplifiedEquation<Kernel, PC, PG>>(); 
            equation->setInputEquation(tg);
            equation->assign(pc);            
            return equation;
        }

        auto equation = std::make_shared<UnaryConstraintComplexEquation<Kernel, PC, PG>>();
        equation->setInputEquation(tg);            
        equation->assign(pc);
        return equation;
    };
    
    virtual std::pair<Equation, FlowNode>
    buildEquationNode(Equation g, symbolic::Constraint* c,
                      shedule::FlowGraph& flowgraph) const override {
       
        auto  tg = std::static_pointer_cast<numeric::Equation<Kernel, PG>>(g);
        auto& pc = static_cast<symbolic::TypeConstraint<PC>*>(c)->getPrimitive();
        
        if(tg->getComplexity() != Complexity::Complex) {
            auto equation = std::make_shared<UnaryConstraintSimplifiedEquation<Kernel, PC, PG>>(); 
            equation->setInputEquation(tg);            
            equation->assign(pc);
            return std::make_pair(equation, flowgraph.newActionNode([=](const shedule::FlowGraph::ContinueMessage& m){
                equation->calculate();
            }));
        }

        auto equation = std::make_shared<UnaryConstraintComplexEquation<Kernel, PC, PG>>();
        equation->setInputEquation(tg);     
        equation->assign(pc);
        return std::make_pair(equation, flowgraph.newActionNode([=](const shedule::FlowGraph::ContinueMessage& m){
            equation->calculate();
        }));
    };
};

}//numeric


//Provide a few default constraints.
//static is needed to restrain the scope of the objects to the current compilation unit. Without it
//every compiled file including this header would define these as global and the linker would find
//multiple definitions of the same objects

struct Distance : public dcm::constraint::Constraint<Distance, 2, double, SolutionSpaces> {
    using Constraint::operator=;
    Distance(){};
    Distance(const double& i, SolutionSpaces s) : Constraint(i,s) {};
    
    double&        distance() {return at<0>();};
    SolutionSpaces solutionSpace() {return at<1>();};
};

struct Orientation : public dcm::constraint::Constraint<Orientation, 2, Orientations> {
    using Constraint::operator=;
    Orientation() {};
    Orientation(const Orientations& i) : Constraint(i) {};
    
    Orientations& orientation() {return at<0>();};
};

struct Angle : public dcm::constraint::Constraint<Angle, 2, double> {
    using Constraint::operator=;
    Angle() {};
    Angle(const double& i) : Constraint(i) {}; 
    
    double& angle() {return at<0>();};
};

struct Fix : public dcm::constraint::Constraint<Fix, 1, Fixables> {
    using Constraint::operator=;
    Fix() {};
    Fix(const Fixables& f) : Constraint(f) {}; 
    
    Fixables& fixed() {return at<0>();};
};

static Distance         distance(0, SolutionSpaces::Bidirectional);
static Orientation      orientation(Orientations::Parallel);
static Angle            angle(0);
static Fix              fix(Fixables::pointX);

}//dcm

#endif //DCM_CONSTRAINT_H





