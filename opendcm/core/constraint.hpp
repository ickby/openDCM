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
#include "geometry.hpp"
#include "defines.hpp"

namespace mpl = boost::mpl;
namespace fusion = boost::fusion;

namespace dcm { 
    
//the possible directions
enum class Orientations { Parallel, Equal, Opposite, Perpendicular };

//the possible solution spaces
enum class SolutionSpaces {Bidirectional, Positiv_directional, Negative_directional};

/**
 * @brief Constraint primitives handling
 *
 * Equivalent to the geometric primitives there is a definition for constraint primitives. A constraint 
 * primitive is a constraint type like distace, primitive in this context means that it holds only its type
 * and all needed information to fully define it, called option. For a distance constraint this is for 
 * example the nuemeric distance value. A primitive is used to store and access the type and its data in a
 * compfortable and user friendly way.
 * There are a few requirements for primitive constraints. First they must be assignalble and copyconstructible
 * while preserving their option values. Second all options need to be stored in a  fusion vector called
 * m_storage. Third it must be possible to provide options via operator() for one or many options at once 
 * and furthermore through assigning with operator=, also for single options or multiple ones via initializer
 * lists. The last requirement regards the default values of the options. It is important to have the possibility
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
 *  void setDefault(); //requirement 4;
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

    /**
     * @brief Herlper class to create primitive constraints
     * 
     * The primitive constraint type concept formulates a few requirements on the used types. To achieve a full
     * compatibility with the concept a certain boilerplate is needed. To circumvent this boilerplate this class
     * is given. It accepts the option types for the constraint as template arguments and then provides
     * all needed storages and assignment operators. 
     * This only works for option types which are default constructible. The only second restriction is that the
     * option types need to distuinguishable. To use a option type twice is not allowd, for example two times int.
     * Otherwise the option assignement of single options would be ambigious. 
     * The default value requirement is achieved through the constructor which acceppts the default values. If 
     * the standart constructor is used then the default constructed values of the option types are used as 
     * default values. Note that you must provide the constructors also in your derived class.
     * An example cosntraint can look like this:
     * @code
     * struct TetstConstraint : public constraint::Constraint<int, char> {
     *    TestConstraint(const int& i, const char& c) : Constraint(i,c) {};
     * };
     * TestConstraint test(1,'a');
     * @endcode
     * 
     * \param OptionTypes a variadic sequence of copy constructable types which describe the stored options
     */    
template<typename ...OptionTypes>
struct Constraint {

    typedef typename fusion::vector<OptionTypes...> Options;
    typedef typename mpl::at_c<Options, 0>::type    PimaryOptionType;

    Constraint() {};
    
    Constraint(const OptionTypes&... defaults) : m_defaults(defaults...) {
        m_storage = m_defaults;
    };
    
    //Copy assign option. We provide this to allow the automatic copy constructor to be generated. 
    //This gives 2 operators for the price of implementing one, for this base and all derived classes.
    //As we have only one parameter in this function it would be ambigious with the single option 
    //operator(), hence we need to activate/deactivate both according to the provided type.
    template<typename T>
    typename boost::disable_if<mpl::contains<Options, T>, Constraint&>::type operator()(const T& c) {
        m_storage = c.m_storage;
        return *this;
    };

    //Set a single option. As we have only one parameter in this function it would be ambigious with the 
    //copy operator(), hence we need to activate/deactivate both according to the provided type.
    template<typename T>
    typename boost::enable_if<mpl::contains<Options, T>, Constraint&>::type operator()(const T& val) {
        BOOST_MPL_ASSERT((mpl::contains<Options, T>));
        *fusion::find<T>(m_storage) = val;
        return *this;
    };
    
    //set multiple options at once.
    Constraint& operator()(const OptionTypes&... val) {
        m_storage = Options(val...);
        return *this;
    };

    //Assign option. Disable it to avoid confusion with the copy assignement operator
    template<typename T>
    typename boost::enable_if<mpl::contains<Options, T>, Constraint&>::type operator=(const T& val) {
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
    
    //access a option by index 
    template<int idx>
    typename fusion::result_of::at_c<Options, idx>::type getOption() {
        return fusion::at_c<idx>(m_storage);
    };
       
protected:
    Options       m_storage;
    const Options m_defaults;
};

}//constraint
    
    
    
namespace symbolic {
    
struct Constraint {
    
    int type;
};

template<typename PrimitiveConstraint>
struct TypeConstraint : public Constraint {
    
    PrimitiveConstraint& getPrimitiveConstraint() {
        return m_constraint;
    };
    
    void setPrimitiveConstraint(const PrimitiveConstraint& c) {
        //type = id;
        m_constraint = c;
    };
    
    void setConstraintID(int id) {
        type = id;
    }
    
protected:   
    PrimitiveConstraint m_constraint;
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
    
}//symbolic    
    
namespace numeric {
   
/**
 * @brief Base class to unify derivation of parent classes
 * 
 * To ease the specification of the parent classes for specilized numeric constraint classes as well 
 * as simplifying the needed types this class is given. It encapsulates all needed inheritance and 
 * typedefs. Note that the equation exposes a single scalar as result. This is due to the fact, that 
 * a constraint equation is always a error function for the numeric solver.
 */
template<typename Kernel, typename PC, template<class> class PG1, template<class> class PG2>
struct ConstraintBase : public BinaryEquation<Kernel, PG1<Kernel>, PG2<Kernel>, typename Kernel::Scalar>, 
                        public PC  {
    
        typedef BinaryEquation<Kernel, PG1<Kernel>, PG2<Kernel>, typename Kernel::Scalar> Equation;
        
        typedef typename Kernel::Scalar                  Scalar;
        typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Vector;
        typedef PG1<Kernel>                              Geometry1;
        typedef Geometry1                                Derivative1;
        typedef PG2<Kernel>                              Geometry2;
        typedef Geometry2                                Derivative2;
};
    
/**
 * @brief Class for numeric evaluation of primitive constraints
 * 
 * Primitive constraints hold only their constraint type and the options to fully define their behaviour. 
 * As numeric evaluation of constraints may be needed in the solving procesm this information is not enough.
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
 * struct Constraint<Kernel, dcm::Distance, dcm::Point3, dcm::Point3> : public ConstraintBase<Kernel, dcm::Distance, TPoint3, TPoint3>  {
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
template<typename Kernel, typename PC, template<class> class PG1, template<class> class PG2>
struct Constraint : public ConstraintBase<Kernel, PC, PG1, PG2> {
           
        typedef ConstraintBase<Kernel, PC, PG1, PG2>     Inherited;
        typedef typename Kernel::Scalar                  Scalar;
        typedef typename Inherited::Vector               Vector;
        typedef typename Inherited::Geometry1            Geometry1;
        typedef typename Inherited::Derivative1          Derivative1;
        typedef typename Inherited::Geometry2            Geometry2;
        typedef typename Inherited::Derivative2          Derivative2;
        
        typedef PC  InheritedConstraint;
        typedef BinaryEquation<Kernel, Geometry1, Geometry2, Scalar> InheritedEquation;
        
        Constraint() {
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
template<typename Kernel, typename PC, template<class> class PG1, template<class> class PG2>
struct ConstraintEquationBase : public numeric::Constraint<Kernel, PC, PG1, PG2> {
   
    typedef Constraint<Kernel, PC, PG1, PG2>            Inherited;
    typedef VectorEntry<Kernel>                         Residual;
    typedef MatrixEntry<Kernel>                         Derivative;
    
    //type to hold geometric derivative together with the correct position for the jacobi entry
    typedef std::pair<typename numeric::Equation<Kernel, PG1<Kernel>>::Derivative*, Derivative> Derivative1Pack;
    typedef std::pair<typename numeric::Equation<Kernel, PG2<Kernel>>::Derivative*, Derivative> Derivative2Pack;
     
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
            g1_derivatives.push_back({&der.first, sys.mapJacobi(residual.Index, der.second.Index)});
    
        for(auto& der : Inherited::secondInputEquation()->derivatives())
            g2_derivatives.push_back({&der.first, sys.mapJacobi(residual.Index, der.second.Index)});
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
#ifdef DCM_DEBUG
        dcm_assert(result1.rows() == g1_derivatives.size());
#endif
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
#ifdef DCM_DEBUG
        dcm_assert(result2.rows() == g2_derivatives.size());
#endif
        int i = 0;
        for(Derivative2Pack& der : g2_derivatives) {
            *(der.second.Value) = result2(i);
            ++i;
        }
    };
    
    void secondAsComplex() {
        
        for(Derivative2Pack& der : g2_derivatives) 
            *(der.second.Value) = Inherited::calculateGradientFirst(Inherited::firstInput(),
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
 * @brief Error function evaluation for two simple inputs
 * 
 * If both input equations are simple, meaning they are not an InputEquation themself, one can use 
 * the simplified functions for both inputs.  This class provides the needed claculation functionality
 * for the equation.
 */
template<typename Kernel, typename PC, template<class> class PG1, template<class> class PG2>
struct ConstraintSimplifiedEquation : public ConstraintEquationBase<Kernel, PC, PG1, PG2> {
    
    typedef ConstraintEquationBase<Kernel, PC, PG1, PG2> Inherited;
    
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
template<typename Kernel, typename PC, template<class> class PG1, template<class> class PG2>
struct ConstraintComplexEquation : ConstraintEquationBase<Kernel, PC, PG1, PG2> {
  
    typedef ConstraintEquationBase<Kernel, PC, PG1, PG2> Inherited;
    
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
template<typename Kernel, typename PC, template<class> class PG1, template<class> class PG2>
struct ConstraintSimplifiedComplexEquation : ConstraintEquationBase<Kernel, PC, PG1, PG2> {
    
    typedef ConstraintEquationBase<Kernel, PC, PG1, PG2> Inherited;
    
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
template<typename Kernel, typename PC, template<class> class PG1, template<class> class PG2>
struct ConstraintComplexSimplifiedEquation : ConstraintEquationBase<Kernel, PC, PG1, PG2> {
  
    typedef ConstraintEquationBase<Kernel, PC, PG1, PG2> Inherited;
   
    CALCULATE() {
        *Inherited::residual.Value = Inherited::calculateError(Inherited::firstInput(), Inherited::secondInput());
        Inherited::firstAsComplex();
        Inherited::secondAsSimplified();
    };
};

template<typename Kernel>
struct ConstraintEquationGenerator {
    
     typedef shedule::FlowGraph::Node  FlowNode;
     typedef std::shared_ptr<numeric::Calculatable<Kernel>>  Equation;
    
     virtual Equation buildEquation(Equation g1, 
                                    Equation g2, 
                                    const symbolic::Constraint* c) const = 0;
                                              
    virtual std::pair<Equation, FlowNode>
    buildEquationNode(Equation g1, 
                      Equation g2, 
                      const symbolic::Constraint* c,
                      shedule::FlowGraph& flowgraph) const = 0;
                                              
};

template<typename Kernel, typename PC, template<class> class PG1, template<class> class PG2>
struct TypedConstraintEquationGenerator : public ConstraintEquationGenerator<Kernel> {

    typedef typename ConstraintEquationGenerator<Kernel>::Equation Equation;
    typedef typename ConstraintEquationGenerator<Kernel>::FlowNode FlowNode;
    
    virtual Equation buildEquation(Equation g1, 
                                   Equation g2, 
                                   const symbolic::Constraint* c) override {
        
        auto tg1 = static_cast<numeric::Equation<Kernel, PG1<Kernel>>*>(g1);
        auto tg2 = static_cast<numeric::Equation<Kernel, PG2<Kernel>>*>(g2);
        
        if(tg1->complexity != Complexity::Complex && tg2->complexity != Complexity::Complex) {
            auto equation = std::make_shared<ConstraintSimplifiedEquation<Kernel, PC, PG1, PG2>>(); 
            equation->setInputEquations(tg1, tg2);            
            return equation;
        }
        else if(tg1->complexity != Complexity::Complex && tg2->complexity != Complexity::Complex) {
            auto equation = std::make_shared<ConstraintSimplifiedComplexEquation<Kernel, PC, PG1, PG2>>();
            equation->setInputEquations(tg1, tg2);            
            return equation;
        }
        else if(tg1->complexity != Complexity::Complex && tg2->complexity != Complexity::Complex) {
            auto equation = std::make_shared<ConstraintComplexSimplifiedEquation<Kernel, PC, PG1, PG2>>();
            equation->setInputEquations(tg1, tg2);            
            return equation;
        }
        else if(tg1->complexity != Complexity::Complex && tg2->complexity != Complexity::Complex) {
            auto equation = std::make_shared<ConstraintComplexEquation<Kernel, PC, PG1, PG2>>();
            equation->setInputEquations(tg1, tg2);            
            return equation;
        }
        dcm_assert(false);
        return Equation();
    };
    
    virtual std::pair<Equation, FlowNode>
    buildEquationNode(Equation g1, Equation g2, const symbolic::Constraint* c,
                      shedule::FlowGraph& flowgraph) const override {
       
        auto  tg1 = static_cast<numeric::Equation<Kernel, PG1<Kernel>>*>(g1);
        auto  tg2 = static_cast<numeric::Equation<Kernel, PG2<Kernel>>*>(g2);
        auto& pc  = static_cast<symbolic::TypeConstraint<PC>*>(c)->getPrimitiveConstraint();
        
        if(tg1->complexity != Complexity::Complex && tg2->complexity != Complexity::Complex) {
            auto equation = std::make_shared<ConstraintSimplifiedEquation<Kernel, PC, PG1, PG2>>(); 
            equation->setInputEquations(tg1, tg2);            
            equation->output() = pc;
            return std::make_pair(equation, flowgraph.newActionNode([=](const shedule::FlowGraph::ContinueMessage& m){
                equation->calculate();
            }));
        }
        else if(tg1->complexity != Complexity::Complex && tg2->complexity != Complexity::Complex) {
            auto equation = std::make_shared<ConstraintSimplifiedComplexEquation<Kernel, PC, PG1, PG2>>();
            equation->setInputEquations(tg1, tg2);     
            equation->output() = pc;
            return std::make_pair(equation, flowgraph.newActionNode([=](const shedule::FlowGraph::ContinueMessage& m){
                equation->calculate();
            }));
        }
        else if(tg1->complexity != Complexity::Complex && tg2->complexity != Complexity::Complex) {
            auto equation = std::make_shared<ConstraintComplexSimplifiedEquation<Kernel, PC, PG1, PG2>>();
            equation->setInputEquations(tg1, tg2);  
            equation->output() = pc;
            return std::make_pair(equation, flowgraph.newActionNode([=](const shedule::FlowGraph::ContinueMessage& m){
                equation->calculate();
            }));
        }
        else if(tg1->complexity != Complexity::Complex && tg2->complexity != Complexity::Complex) {
            auto equation = std::make_shared<ConstraintComplexEquation<Kernel, PC, PG1, PG2>>();
            equation->setInputEquations(tg1, tg2);  
            equation->output() = pc;
            return std::make_pair(equation, flowgraph.newActionNode([=](const shedule::FlowGraph::ContinueMessage& m){
                equation->calculate();
            }));
        }
        dcm_assert(false);
        return Equation();
    };
};

}//numeric


//Provide a few default constraints.
//static is needed to restrain the scope of the objects to the current compilation unit. Without it
//every compiled file including this header would define these as global and the linker would find
//multiple definitions of the same objects

struct Distance : public dcm::constraint::Constraint<double, SolutionSpaces> {
    using Constraint::operator=;
    Distance(){};
    Distance(const double& i, SolutionSpaces s) : Constraint(i,s) {};
    
    double&        distance() {return fusion::at_c<0>(m_storage);};
    SolutionSpaces solutionSpace() {return fusion::at_c<1>(m_storage);};
};

struct Orientation : public dcm::constraint::Constraint<Orientations> {
    using Constraint::operator=;
    Orientation(){};
    Orientation(const Orientations& i) : Constraint(i) {};
    
    Orientations& orientation() {return fusion::at_c<0>(m_storage);};
};

struct Angle : public dcm::constraint::Constraint<double> {
    using Constraint::operator=;
    Angle(){};
    Angle(const double& i) : Constraint(i) {};
    
    double& angle() {return fusion::at_c<0>(m_storage);};
};

static Distance         distance(0, SolutionSpaces::Bidirectional);
static Orientation      orientation(Orientations::Parallel);
static Angle            angle(0);

}//dcm

#endif //DCM_CONSTRAINT_H





