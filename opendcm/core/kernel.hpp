/*
    openDCM, dimensional constraint manager
    Copyright (C) 2012  Stefan Troeger <stefantroeger@gmx.net>

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

#ifndef DCM_KERNEL_H
#define DCM_KERNEL_H

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <boost/graph/graph_concepts.hpp>

#include "transformation.hpp"
#include "logging.hpp"
#include "scheduler.hpp"

namespace dcm {
namespace numeric {
    
template<typename Kernel> 
struct LinearSystem;

  
//every kernel needs to be derived from this class
struct KernelBase {};
    
template< typename Kernel >
struct VectorEntry {
    
    typedef typename Kernel::Scalar Scalar;
    
    int     Index;
    Scalar* Value;
    
    bool operator==(const VectorEntry& e) {
        //comparing the value addresses is enough, as if tey point to the same memory they
        //must be the same system entry
        return (Value == Value);
    };
    
    operator Scalar&() {return *Value;};
};

template< typename Kernel >
struct MatrixEntry {
    
    int                         Row;
    int                         Column;
    typename Kernel::Scalar*    Value;
    
    bool operator==(const MatrixEntry& e) {
        //comparing the value addresses is enough, as if tey point to the same memory they
        //must be the same system entry
        return (Value == Value);
    };
};

template<typename Kernel> 
struct LinearSystem {
    
    typedef typename Kernel::Scalar                  Scalar;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> VectorX;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic,
                Eigen::Dynamic>                      MatrixX;
    
    LinearSystem(int p, int e) : m_parameterCount(p), m_equationCount(e),
            m_jacobi(e, p), m_parameters(p),  m_residuals(e) {};
    
    
    VectorEntry<Kernel> mapParameter() {
        Scalar* s = &m_parameters(++m_parameterOffset);
        return {m_parameterOffset, s};
    };                
               
    std::vector<VectorEntry<Kernel>> mapParameter(Scalar*& s) {
        s = &m_parameters(++m_parameterOffset);
        std::vector<VectorEntry<Kernel>> result(1);
        result[0] = {m_parameterOffset, s};
        return result;
    };

    template<typename Derived>
    std::vector<VectorEntry<Kernel>> mapParameter(Eigen::Map<Derived>& map) {
        new(&map) Eigen::Map<Derived>(&m_parameters(++m_parameterOffset));
        
        std::vector<VectorEntry<Kernel>> result(map.rows());
        for (int i = 0; i < map.rows(); ++i)
            result[i] = {m_parameterOffset + i, &m_parameters(m_parameterOffset + i)};

        //minus one as we increase the parameter offset already in this function
        m_parameterOffset+= (map.rows()-1);
                
        return result;
    };
    
    VectorEntry<Kernel> mapResidual(Scalar*& s) {
        s = &m_residuals(++m_residualOffset);
        return {m_residualOffset, s};
    };
    
    VectorEntry<Kernel> mapResidual() {
        Scalar* s = &m_residuals(++m_residualOffset);
        return {m_residualOffset, s};
    };
    
    MatrixEntry<Kernel> mapJacobi(int row, int col, Scalar*& s) {
        s = &m_jacobi(row, col);
        return {row, col, s};
    };  
    
    MatrixEntry<Kernel> mapJacobi(int row, int col) {
        Scalar* s = &m_jacobi(row, col);
        return {row, col, s};
    };   
    
    void setupJacobi() {
        m_jacobi.setZero();
    };
    
    Scalar& jacobiAt(int row, int col) {
       return m_jacobi(row, col);  
    };
    
    //access the vectors and matrices
    VectorX& parameter() {return m_parameters;};
    VectorX& residuals() {return m_residuals;};
    MatrixX& jacobi()    {return m_jacobi;};    
    
private:
    int m_parameterCount, m_equationCount;
    int m_parameterOffset = -1, m_residualOffset  = -1;
    VectorX m_parameters;
    VectorX m_residuals;
    MatrixX m_jacobi;
};


/**
 * @brief Base for all things that need to be calculated when evaluating a Linear System
 * 
 */
template<typename Kernel>
struct Calculatable : shedule::Executable {
    
    typedef Kernel KernelType;
    
    /**
     * @brief Initialization to enable internal setup
     * 
     * This function must be called before any access to the calculatable is made. It ensures that all inernals are
     * valid, all parameters are initialised. If the calculatabe is accessed before calling init the behaviour can be
     * undefined. In debug mode an assert will be called in this situation, but in release no warning will occure. 
     * The geometry can only be initialized with a \ref LinearSystem as the parameters of the equation are maped
     * into this lienar system. It is therefore highly important that the given \ref LinearSystem is used for
     * all calculations involing this class.
     * 
     * @param sys LinearSystem the geometry is initalized with
     * @return void
     */
    virtual void init(LinearSystem<Kernel>& sys) {};
    
    /**
     * @brief Execution of the calculation
     * 
     * This function executes the calculation. It is the derived classes responsibility to override this 
     * method and provide suitable behaviour
     * 
     * The caller is responsible for ensuring that the equation has been initialized before call this function.
     * @return void
     */
    virtual void execute() {};  
    
    /**
     * @brief Number of free parameters this equation needs
     * 
     * Return the amount of parameters this equation will need in the solving process. The amount is later
     * mapped into the linear system on \ref init. It is possible to access this function before initialisation.
     * @return void
     */
    unsigned int newParameterCount() {return m_parameterCount;};
    
protected:
    int m_parameterCount = 0; //how many parameters are added by this equation?

};
    
//the standart solverS
template<typename Kernel>
struct Dogleg {

#ifdef DCM_USE_LOGGING
    dcm_logger log;
#endif

    typedef typename Kernel::Scalar Scalar;
    Scalar tolg, tolx, delta, nu, g_inf, fx_inf, err, time;
    Kernel* m_kernel;
    int iter, stop, reduce, unused, counter;
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> h_dl, F_old, g;
    Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> J_old;

    Dogleg(Kernel* k);
    Dogleg();
/*
    void setKernel(Kernel* k);

    template <typename Derived, typename Derived2, typename Derived3, typename Derived4>
    void calculateStep(const Eigen::MatrixBase<Derived>& g, const Eigen::MatrixBase<Derived3>& jacobi,
                       const Eigen::MatrixBase<Derived4>& residual, Eigen::MatrixBase<Derived2>& h_dl,
                       const double delta);

    int solve(typename Kernel::MappedEquationSystem& sys);

    template<typename Functor>
    int solve(typename Kernel::MappedEquationSystem& sys, Functor& rescale);*/
};

struct DummyKernel : public numeric::KernelBase {

    typedef int Scalar;
};

};

template<typename NumericType, template<class> class Nonlinear = numeric::Dogleg>
struct Eigen3Kernel : public numeric::KernelBase {

    //the number type we use throughout the system
    typedef NumericType   Scalar;


private:
    Nonlinear< Eigen3Kernel<Scalar, Nonlinear> > m_solver;

};

}//dcm

#ifndef DCM_EXTERNAL_CORE
//#include "imp/kernel_imp.hpp"
#endif

#endif //GCM_KERNEL_H






