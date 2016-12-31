/*
    openDCM, dimensional constraint manager
    Copyright (C) 2016  Stefan Troeger <stefantroeger@gmx.net>

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

#ifndef DCM_EQUATIONS_H
#define DCM_EQUATIONS_H

#include <functional>
#include <memory>

#include <boost/mpl/if.hpp>
#include <boost/type_traits.hpp>

#include "defines.hpp"
#include "kernel.hpp"

namespace dcm {

namespace detail {
    
    

/**
    * @brief Allow inheritance from POD types
    * 
    * Some template classes, as equations, are derived from user provided classes. Sometimes it is 
    * possible that those classes are c++ pod types and therefore one can not use them for inheritance.
    * This helper class wraps the pod and makes it accessible by inheritance. It therefore allows 
    * assigning, copying and casting.
    * 
    * \tparam Pod the c++ pod type that should be made inheritable
    */
template<typename Pod>
struct PodBase {
  
    PodBase() {};
    PodBase(const Pod& val) : m_pod(val) {};
    
    Pod& value() {
        return m_pod;
    }
    
    operator Pod&() {
        return m_pod;
    }
    
    void operator=(const Pod& val) {
        m_pod = val;
    }
    
private:
    Pod m_pod;
};

} //detail

namespace numeric {
    
namespace mpl = boost::mpl;

//Fixed:   No free parameters, fixed output value
//Simple:  Every value of the output is an free parameter, equation depends on no additional ones
//Complex: Output is calculated from free parameters of the own or other equations
enum class Complexity { Fixed, Simple, Complex };


//This macro is intendet to ease the implementation of equations. They need to override execute() to allow 
//polymorphic calculation, but also implement calculate() for exeution without virtual calls. Futhermore 
//using them as functors in heplpful. This macro provides the boilerplate needed for all of this.
#define CALCULATE() \
    virtual void execute() {\
        calculate();\
    }; \
    void operator()() {\
        calculate(); \
    }; \
    void calculate()
               
/**
    * @brief Base class for numeric equations
    * 
    * This class provides the common infrastructure for calculations of equations. An dcm equation 
    * does not mean the same as the mathematical understanding of f(x), it is the whole 
    * expression y = f(x). This means it exposes the result y and provides the calculation f(x). 
    * Furthermore it is responsible to calculate the derivative y' = f'(x) too. 
    * 
    * The class provides the storage and access for the free parameters it needs as well as the 
    * derivatives of the result from all relevant parameters. As both is not unique but implementation 
    * dependend no values are assigned in this base class. 
    * 
    * \ref Equation doe inherit the \ref Output template parameter. This means it can directly be used
    * as result y as well as equation f(x).
    * 
    * \note This class can be used as a fixed value equation when given the value to the constructor. 
    *       It does than not have any free parameters or derivatives but can be usefull as input for
    *       other equations.
    * 
    * \tparam Kernel The mathematical kernel to be used 
    * \tparam Output The result type y that the equation calculates
    */
template<typename Kernel, typename Output>
struct Equation : public mpl::if_<boost::is_pod<Output>, detail::PodBase<Output>, Output>::type, 
                  public Calculatable<Kernel>,
                  public std::enable_shared_from_this<Equation<Kernel, Output>> {   
    
    typedef typename mpl::if_<boost::is_pod<Output>, detail::PodBase<Output>, Output>::type Base;
    using Base::operator=;
    
    typedef Kernel                                              KernelType;
    typedef Output                                              OutputType;
    typedef VectorEntry<Kernel>                                 Parameter;
    typedef typename std::vector<Parameter>::iterator           ParameterIterator;
    typedef Output                                              Derivative;
    typedef std::pair<Output, Parameter>                        DerivativePack;
    typedef typename std::vector<DerivativePack>::iterator      DerivativePackIterator;
    
    //set values in case this is a fixed equation
    Equation() {};
    Equation(const Output& val) : Base(val) {};
    
    /**
     * @brief Gives the result
     * 
     * The equation derives from the result type and can therefore be directly used as output. However,
     * it sometimes may be nessecary to cast the equation directly to the result type. This is done 
     * by this function for direct access.
     * @return Output& the result of the equations calculation
     */
    Output& output() {return *static_cast<Base*>(this);}
   
    /**
     * @brief Cast equation to result type
     * The equation derives from the result type and can therefore be directly used as output. However,
     * it sometimes may be nessecary to cast the equation directly to the result type.
     * \see output()
     * @return const Output&
     */    
    operator const Output&() const {return output();} 
    
    Output& operator=(const Output& in) {output() = in;};
    
    /**
     * @brief Access all initialized free parameters
     * 
     * Allows to access all the free parameters this equation adds to the overall system. Note that 
     * Only the added parameters are returned here. The equation may depend on other parameters as 
     * well, for example from input equations, but they are not listed here. Hence it can also 
     * happen that there are more derivatives available than parameters.
     * 
     * @return std::vector< Parameter >& vector of all mapped free parameters
     */
    std::vector< Parameter >&      parameters() {
        return m_parameters;
    };    
    
    /**
     * @brief Access all initialized derivatives of this equation
     * 
     * This equation may depend on multiple parameters, the ones we allocated ourself (see \ref parameter) as 
     * well as allocate by other equations we depend on. The returned values are only valid after initialisation.
     * The returned vector allows to access the parameter for which the derivative is calculated 
     * as well as the derivative itself, which is given by the same type as the result itself is. 
     * 
     * @return std::vector< Parameter >& vector of all mapped free parameters
     */
    std::vector< DerivativePack >& derivatives() {
        return m_derivatives;
    };
    
    /**
     * @brief Returns the complexity of this equation
     * 
     * Returns a hint on how this equation is calculated:
     * Fixed:   No free parameters, fixed output value
     * Simple:  Every value of the output is an free parameter, equation depends on no additional ones
     * Complex: Output is calculated from free parameters of the own or other equations
     * 
     * @return dcm::numeric::Complexity Complexity of the equation
     */
    Complexity getComplexity() {return m_complexity;};
    
#ifdef DCM_DEBUG
    bool isInitialized() {
        return m_init;
    };
#endif
    
protected:   
    //storage of derivatives for faster calculation
    std::vector< Parameter >            m_parameters;
    std::vector< DerivativePack >       m_derivatives; 
    Complexity                          m_complexity = Complexity::Fixed;
#ifdef DCM_DEBUG
    bool m_init = false;
#endif
    
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

/**
 * @brief Handles ownership of equations with inputs 
 * 
 * For equations with inputs other than the free parameters, e.g. other equations, the question of ownership arises.
 * There may be the possibility that the input equation is initialised and calculated by someone else, than this 
 * does not need to be redone by the equation using it. However, sometimes it may happen that the input equation 
 * is owned by the equation itself, hence its init and calculate methods will not be called by annyone else. In
 * this case the equation is responsible for handling this. \ref this class provides the interface needed to 
 * tell the equation if it owns its input equations and if it is responsible for handling it.
 */
template<typename Kernel, typename Output>
struct InputEquation : public Equation<Kernel, Output> {
    
    InputEquation() {
        Equation<Kernel, Output>::m_complexity = Complexity::Complex;
    };
    
    /**
     * @brief Returns if equation owns the input equations
     * 
     * If the input equations are nether initialized or executed from someone and this equation 
     * is responsible for doing so, this function will return true.
     * 
     * @return bool True if it owns the input equation, false otherwise
     */
    bool hasInputOwnership() {return m_ownership;}
    
    /**
     * @brief Set ownership of input equations
     * 
     * Allows to specify if this equation owns the input equations and should be responsible for 
     * initializing and calculating them. 
     * 
     * \param val True if ownership should be assigned to this equation
     */
    void takeInputOwnership(bool val) {m_ownership = val;};
    
protected:
    bool m_ownership = false;
};



/**
 * @brief Equation depending on one inut equation
 * 
 * This equation extends the base class with all functionality needed for handling a single input 
 * equation. Therefore it is not only defined by the result type y it outputs, like the base equation, 
 * but also by the input type j it needs for its calculation: y = f(x,j); Therefore only equations with
 * an output type equal this unary equations input type can be used as input. This also means that the 
 * derivative of y is dependend on the free parameter x as well as all derivatives from input j. Hence
 * the derivatives vector will be larger than the parameters one.
 * 
 * This class provides functions for easy handling of external owned equations as input as well as 
 * internall owned ones. 
 * 
 * \tparam Kernel The mathematical kernel in use 
 * \tparam Input The input type needed for the calculation
 * \tparam Output The result type provided by this equation
 */
template<typename Kernel, typename Input, typename Output>
struct UnaryEquation : public InputEquation<Kernel, Output> {
    
    typedef Input                               InputType;
    typedef Equation<Kernel, Input>             InputEqn;
    typedef InputEquation<Kernel, Output>       Base;
    
    UnaryEquation() {};
    UnaryEquation(std::shared_ptr<InputEqn> in) : m_input(in) {};
    
    /**
     * @brief Access the input value for this unary equation
     * 
     * Returns the input value for this equation, which is basially the output of the equation used 
     * as input. 
     * \note This function can only be called after the input equation has been set.
     * 
     * \return Input& The value of the input used for calculations or the output
     */
    Input&  input() {
        dcm_assert(m_input);
        return m_input->output();       
    };
    
    /**
     * @brief Access the input equation for this unary equation
     * 
     * Returns the Equation used as input for this equation. This is the one previously set.
     * \note This function returns an NULL equation if called before the equation was set
     * @return std::shared_ptr< InputEqn >
     */    
    std::shared_ptr<InputEqn>   inputEquation() {return m_input;};
    
    /**
     * @brief Set the input equation 
     * 
     * Allows to specify the input equation used by this unary equation. If this function is used it 
     * is assumed that the ownership of the \ref eqn lies outside of \ref this. This means 
     * \ref hasInputOwnership will return false and the inputs init and calculate functions won't be 
     * called by this equation.
     * @return void
     */
    void  setInputEquation(std::shared_ptr<InputEqn> eqn) {
        m_input = eqn;
        Base::takeInputOwnership(false);
    }
    
    
    /**
     * @brief Append unary equation to this equation
     * 
     * This functions chains two unary equations. The given parameter is appended to \a this, 
     * which makes \a this the input equation of the given equaiton \ref ptr. This function does also 
     * transfer ownership for \a this equation to \ref ptr. 
     * The equation to append must have the same input type as \a this has as output type. 
     * \param ptr The UnaryEquation to append to \a this
     * \return The Equation providing the new output
     */
    template<typename NewOutput> 
    std::shared_ptr<Equation<Kernel, NewOutput>> 
    append(std::shared_ptr<UnaryEquation<Kernel, Output, NewOutput>> ptr) {
        
        ptr->setInputEquation(Base::shared_from_this());
        ptr->takeInputOwnership(true);
        return ptr;
    };
    
    /**
     * @brief Prepend an equation to \a this
     * 
     * This functions chains two equations. The given parameter \ref ptr is prepended to \a this, 
     * which makes it the input equation of \a this. This function does also transfer ownership for
     * \ref ptr equation to \a this. The equation to prepend must have the same output type as 
     * \a this has as input type. 
     * \param ptr The Equation to prepend to \a this
     * \return The Equation providing the new output
     */
    std::shared_ptr<Equation<Kernel, Output>> 
    prepend(std::shared_ptr<Equation<Kernel, Input>> ptr) {
        
        setInputEquation(ptr);
        Base::takeInputOwnership(true);
        return Base::shared_from_this();
    };
    
protected:
    std::shared_ptr<InputEqn> m_input;
};

/**
 * @brief UnaryEquation basen  on expressions
 * 
 * This class is a fully defiend unary expression where all calculations are defined by Expressions. 
 * To calculate the result two expressions are needed: one for the result and one for the derivative.
 * The derivative expression is called as often as many parameter this equation depends on. 
 * 
 * An expression is an arbitrary functor object with a special syntax for calculations. The structure
 * of the expression for the function evaluation is given below as lambda example.
 * \code{.cpp}
 * [](const Input& in, Output& out) {
*      out = 2*in*in;    
*  },        
 * \endcode
 * The derivatives are calculates as following:
 * \code{.cpp}
 * [](const Input& in, const Input& derivative_in, Output& out) {
 *     out = 4*in*derivative_in;
 * }
 * \endcode
 * 
 * There is no init expression supportet, hence an ExpressionUnaryEquation cannot have free parameters.
 * 
 * \note The expressions passed to the constructor are copyed, hence using members by value is a bad 
 *       idea.  
 */
template<typename Kernel, typename Input, typename Output, typename CExp, typename DExp>
struct ExpressionUnaryEquation : public UnaryEquation<Kernel, Input, Output> {
    
    typedef UnaryEquation<Kernel, Input, Output> Base;
    
    ExpressionUnaryEquation(const CExp& c, const DExp& d) : m_cExp(c), m_dExp(d) {}
    
    virtual void init(LinearSystem<Kernel>& k) {
        
        dcm_assert(Base::m_input);
        if(Base::hasInputOwnership())
            Base::m_input->init(k);
               
        //copy over the derivative parameters. As we don't add any new parameters this is sufficient
        Base::m_derivatives.clear();
        for(const auto& param : Base::inputEquation()->parameters()) 
            Base::m_derivatives.push_back(std::make_pair(Output(), param));
    }
    
    CALCULATE() {
        
        dcm_assert(Base::m_input);        
        if(Base::hasInputOwnership())
            Base::m_input->execute(); 
        
        m_cExp(Base::input(), Base::output());
        
        auto& vec = Base::m_derivatives;
        typename std::vector< typename Base::DerivativePack >::iterator it = vec.begin();
        
        for(auto& der : Base::inputEquation()->derivatives()) {
            m_dExp(Base::input(), der.first, it->first);
            ++it;
        }
    };
    
private:
    CExp m_cExp;
    DExp m_dExp;
};

/**
 * @brief Creates a UnaryEquation from expressions
 * 
 * Returns a fully qualified UnaryEquation based on the given parameters. This is achieved by creating
 * a \ref  ExpressionUnaryEquation based on the parameters. Hence the expressions need to have the structure
 * defined by \ref ExpressionUnaryEquation. Normal usage is with lambda expressions. Note that the input 
 * and output types need to be spezified:
 * \code{.cpp}
 * auto d_v = numeric::makeUnaryEquation<K, double, Eigen::Vector3d>  (
 *     [](const double& d, Eigen::Vector3d& v) {
 *         v << d, 2*d, 3*d;    
 *     },        
 *     [](const double& d, const double& dd, Eigen::Vector3d& v) {
 *         v << 1*dd,2*dd,3*dd;
 *     }
 * );
 * \endcode
 */
template<typename Kernel, typename Input, typename Output, typename CExpr, typename DExpr>
std::shared_ptr<UnaryEquation<Kernel, Input, Output>> makeUnaryEquation(const CExpr& cexpr, const DExpr& dexpr) {

    auto ptr = new ExpressionUnaryEquation<Kernel, Input, Output, CExpr, DExpr>(cexpr, dexpr);
    return std::shared_ptr<UnaryEquation<Kernel, Input, Output>>(ptr);
}

/**
 * @brief Creates unary equation dependend on input
 * 
 * This function creates a ExpressionUnaryEquation equal to \ref makeUnaryEquation, but additionally
 * appends it to the given equation. This allows to use this function to "extend" a given equation
 * with expressions.
 */
template<typename Kernel, typename Input, typename Output, typename CExpr, typename DExpr>
std::shared_ptr<UnaryEquation<Kernel, Input, Output>> makeUnaryEquation(std::shared_ptr<Equation<Kernel, Input>> eqn, 
                                                            const CExpr& cexpr, const DExpr& dexpr) {

    auto ptr = makeUnaryEquation<Kernel, Input, Output>(cexpr, dexpr);
    ptr->prepend(eqn);
    return ptr;
}


/**
 * @brief Equation depending on two input equations
 * 
 * This class is an equivalent to the \ref UnaryEquation but with two input equations. It extends the
 * base class with all functionality needed for handling of two input equations. Therefore it is not 
 * only defined by the result type y it outputs, like the base equation,  but also by the input type
 * j and k it needs for its calculation: y = f(x,j, k); Therefore only equations with an output type 
 * equal one of this unary equations input typea can be used as input. This also means that the 
 * derivative of y is dependend on the free parameter x as well as all derivatives from inputa j and 
 * k. Hence the derivatives vector will be larger than the parameters one.
 * 
 * This class provides functions for easy handling of external owned equations as input as well as 
 * internall owned ones. 
 * 
 * \tparam Kernel The mathematical kernel in use 
 * \tparam Input1 The first input type needed for the calculation
 * \tparam Input2 The second input type needed for the calculation
 * \tparam Output The result type provided by this equation
 */
template<typename Kernel, typename Input1, typename Input2, typename Output>
struct BinaryEquation : public InputEquation<Kernel, Output> {
      
    typedef Input1                              Input1Type;
    typedef Input2                              Input2Type;
    typedef Equation<Kernel, Input1>            Input1Eqn;
    typedef Equation<Kernel, Input2>            Input2Eqn;
    typedef InputEquation<Kernel, Output>       Base;
    
    BinaryEquation() {};
    BinaryEquation(std::shared_ptr<Input1Eqn> in1, std::shared_ptr<Input2Eqn> in2) 
        : m_input1(in1), m_input2(in2) {};
    
    /**
     * @brief Access the first input value for this binary equation
     * 
     * Returns the first input value for this equation, which is basially the output of the first 
     * equation used as input. 
     * \note This function can only be called after the first input equation has been set.
     * 
     * \return Input1& The value of the input used for calculations or the output
     */
    Input1&  firstInput() {
        dcm_assert(m_input1);
        return m_input1->output();
    };
    
    /**
     * @brief Access the second input value for this binary equation
     * 
     * Returns the second input value for this equation, which is basially the output of the second 
     * equation used as input. 
     * \note This function can only be called after the second input equation has been set.
     * 
     * \return Input2& The value of the input used for calculations or the output
     */
    Input2& secondInput() {
        dcm_assert(m_input2);
        return m_input2->output();
    };
    
    /**
     * @brief Access the first input equation for this binary equation
     * 
     * Returns the first equation used as input for this equation. This is the one previously set.
     * \note This function returns an NULL equation if called before the equation was set
     * @return std::shared_ptr< Input1Eqn > the first input equation
     */ 
    std::shared_ptr<Input1Eqn>  firstInputEquation() {return m_input1;};
    
    /**
     * @brief Access the second input equation for this binary equation
     * 
     * Returns the second equation used as input for this equation. This is the one previously set.
     * \note This function returns an NULL equation if called before the equation was set
     * @return std::shared_ptr< Input2Eqn > the second input equation
     */ 
    std::shared_ptr<Input2Eqn>  secondInputEquation() {return m_input2;};
    
    /**
     * @brief Set the first input equation 
     * 
     * Allows to specify the first input equation used by this binary equation. If this function is used it 
     * is assumed that the ownership of the \ref eqn lies outside of \ref this. This means 
     * \ref hasInputOwnership will return false and the inputs init and calculate functions won't be 
     * called by this equation.
     */
    void setFirstInputEquation(std::shared_ptr<Input1Eqn> eqn) {m_input1 = eqn;}
    
    /**
     * @brief Set the second input equation 
     * 
     * Allows to specify the second input equation used by this binary equation. If this function is used it 
     * is assumed that the ownership of the \ref eqn lies outside of \ref this. This means 
     * \ref hasInputOwnership will return false and the inputs init and calculate functions won't be 
     * called by this equation.
     */
    void setSecondInputEquation(std::shared_ptr<Input2Eqn> eqn) {m_input2 = eqn;}
    
    /**
     * @brief Set the all input equations
     * 
     * Allows to specify both input equations used by this binary equation. If this function is used it 
     * is assumed that the ownership of the \ref eqnq and \ref eqn2 lies outside of \ref this. This means 
     * \ref hasInputOwnership will return false and the inputs init and calculate functions won't be 
     * called by this equation.
     */
    void setInputEquations(std::shared_ptr<Input1Eqn> eqn1, 
                           std::shared_ptr<Input2Eqn> eqn2) {
                     m_input1 = eqn1; m_input2 = eqn2;}
    
    
    /**
     * @brief Append unary equation to this equation
     * 
     * This functions chains a unary equations to the output of \a this binary equation. The given 
     * parameter \ref ptr is appended to \a this, which makes \a this the input equation of the given 
     * equaiton \ref ptr. This function does also transfer ownership for \a this equation to \ref ptr. 
     * The equation to append must have the same input type as \a this has as output type. 
     * \param ptr The UnaryEquation to append to \a this
     * \return The Equation providing the new output
     */
    template<typename NewOutput> 
    std::shared_ptr<Equation<Kernel, NewOutput>> 
    append(std::shared_ptr<UnaryEquation<Kernel, Output, NewOutput>> ptr) {
        
        ptr->setInputEquation(Base::shared_from_this());
        ptr->takeInputOwnership(true);
        return ptr;
    };
    
    /**
     * @brief Prepend two equations to \a this
     * 
     * This functions chains two equations as inputs to this binary equation. The given parameters 
     * \ref ptr1 and \ref ptr2 are prepended to \a this, which makes them the input equations of \a this.
     * This function does also transfer ownership for \ref ptr1 and \ref ptr2 equations to \a this. 
     * The equations to prepend must have the same output type as \a this has as input types, one each 
     * matching one input type.
     * \param ptr1 The Equation to prepend to the first input of \a this
     * \param ptr2 The Equation to prepend to the second input of \a this
     * \return The Equation providing the new output
     */
    std::shared_ptr<Equation<Kernel, Output>> 
    prepend(std::shared_ptr<Equation<Kernel, Input1>> ptr1, std::shared_ptr<Equation<Kernel, Input2>> ptr2) {
        
        setFirstInputEquation(ptr1);
        setSecondInputEquation(ptr2);
        Base::takeInputOwnership(true);
        return Base::shared_from_this();
    };
    
protected:
    std::shared_ptr<Input1Eqn> m_input1;
    std::shared_ptr<Input2Eqn> m_input2;
};

/**
 * @brief BinaryEquation basen on expressions
 * 
 * This class is a fully defiend binary equation where all calculations are defined by expressions. 
 * To calculate the result three expressions are needed: one for the result and one for each input 
 * derivative. The derivative expressiona are called as often as many parameter this equation depends on. 
 * 
 * An expression is an arbitrary functor object with a special syntax for calculations. The structure
 * of the expression for the function evaluation is given below as lambda example.
 * \code{.cpp}
 * [](const Input& in1, const Input2& in2, Output& out) {
*      out = 2*in1*in2;    
*  },        
 * \endcode
 * The two derivative expressions are calculates as following:
 * \code{.cpp}
 * [](const Input1& in1, const Input2& in2, const Input1& derivative_in1, Output& out) {
 *     out = 2*derivative_in1*in2;
 * }
 * [](const Input1& in1, const Input2& in2, const Input2& derivative_in2, Output& out) {
 *     out = 2*in1*derivative_in2;
 * }
 * \endcode
 * 
 * There is no init expression supportet, hence an ExpressionBinaryEquation cannot have free parameters.
 * 
 * \note The expressions passed to the constructor are copyed, hence using members by value is a bad 
 *       idea.  
 */
template<typename Kernel, typename Input1, typename Input2, typename Output, 
         typename CExp, typename DExp1, typename DExp2>
struct ExpressionBinaryEquation : public BinaryEquation<Kernel, Input1, Input2, Output> {
    
    typedef BinaryEquation<Kernel, Input1, Input2, Output> Base;
    
    ExpressionBinaryEquation(const CExp& c, const DExp1& d1, const DExp2& d2) 
        : m_cExp(c), m_dExp1(d1), m_dExp2(d2) {}
    
    virtual void init(LinearSystem<Kernel>& k) {
        
        dcm_assert(Base::m_input1);
        dcm_assert(Base::m_input2);
        if(Base::hasInputOwnership()) {
            Base::firstInputEquation()->init(k);
            Base::secondInputEquation()->init(k);
        }

        //we add no own parameters
        Base::m_parameters.clear();
        
        //copy over the derivative parameters. As we don't add any new parameters this is sufficient
        Base::m_derivatives.clear();
        for(const auto& param : Base::firstInputEquation()->parameters()) 
            Base::m_derivatives.push_back(std::make_pair(Output(), param));
        
        for(const auto& param : Base::secondInputEquation()->parameters()) 
            Base::m_derivatives.push_back(std::make_pair(Output(), param));
    }
    
    CALCULATE() {
        
        dcm_assert(Base::m_input2);        
        dcm_assert(Base::m_input1); 
        if(Base::hasInputOwnership())  {
            Base::m_input1->execute(); 
            Base::m_input2->execute(); 
        }
        
        m_cExp(Base::firstInput(), Base::secondInput(), Base::output());
        
        auto& vec = Base::m_derivatives;
        typename std::vector< typename Base::DerivativePack >::iterator it = vec.begin();
        
        for(auto& der : Base::firstInputEquation()->derivatives()) {
            m_dExp1(Base::firstInput(), Base::secondInput(), der.first, it->first);
            ++it;
        }
        for(auto& der : Base::secondInputEquation()->derivatives()) {
            m_dExp2(Base::firstInput(), Base::secondInput(), der.first, it->first);
            ++it;
        }
    };
    
private:
    CExp  m_cExp;
    DExp1 m_dExp1;
    DExp2 m_dExp2;
};


/**
 * @brief Creates a BinaryEquation from expressions
 * 
 * Returns a fully qualified BinaryEquation based on the given parameters. This is achieved by creating
 * a \ref  ExpressionBinaryEquation based on the parameters. Hence the expressions need to have the structure
 * defined by \ref ExpressionBinaryEquation. Normal usage is with lambda expressions. Note that both inputs
 * types and the output type need to be spezified:
 * \code{.cpp}
 * auto d_v = numeric::makeUnaryEquation<K, double, double, Eigen::Vector3d>  (
 *     [](const double& d1, const double& d2, Eigen::Vector3d& v) {
 *         v << d1, 2*d1, 3*d2;    
 *     },        
 *     [](const double& d1, const double& d2, const double& dd1, Eigen::Vector3d& v) {
 *         v << 1*dd1,2*dd1,3*0;
 *     },        
 *     [](const double& d1, const double& d2, const double& dd2, Eigen::Vector3d& v) {
 *         v << 1*0,2*0,3*dd2;
 *     }
 * );
 * \endcode
 */
template<typename Kernel, typename Input1, typename Input2, typename Output, 
         typename CExpr, typename DExpr1, typename DExpr2>
std::shared_ptr<BinaryEquation<Kernel, Input1, Input2, Output>> makeBinaryEquation(const CExpr& cexpr, 
                                                                                             const DExpr1& dexpr1,
                                                                                             const DExpr2& dexpr2 ) {

    return std::shared_ptr<BinaryEquation<Kernel, Input1, Input2, Output>>(
            new ExpressionBinaryEquation<Kernel, Input1, Input2, Output, CExpr, DExpr1, DExpr2>(cexpr, dexpr1, dexpr2));
}

/**
 * @brief Creates binary equation dependend on input
 * 
 * This function creates a ExpressionBinaryEquation equal to \ref makeBinaryEquation, but additionally
 * appends itself to the given equations. This allows to use this function to "extend" given equations
 * with expressions.
 */
template<typename Kernel, typename Input1, typename Input2, typename Output, 
         typename CExpr, typename DExpr1, typename DExpr2>
std::shared_ptr<BinaryEquation<Kernel, Input1, Input2, Output>> makeBinaryEquation(
                                                            std::shared_ptr<Equation<Kernel, Input1>> eqn1, 
                                                            std::shared_ptr<Equation<Kernel, Input2>> eqn2, 
                                                            const CExpr& cexpr, const DExpr1& dexpr1,
                                                            const DExpr2& dexpr2) {

    auto ptr = makeBinaryEquation<Kernel, Input1, Input2, Output>(cexpr, dexpr1, dexpr2);
    ptr->prepend(eqn1, eqn2);
    return ptr;
}


} //numeric    
} //dcm

#endif //DCM_EQUATIONS_H


