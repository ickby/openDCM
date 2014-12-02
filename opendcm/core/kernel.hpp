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

namespace dcm {
namespace numeric {
    
//every kernel needs to be derived from this class
struct KernelBase {};
    
template< typename Kernel >
struct SystemEntry {
    
    int                         Index;
    typename Kernel::Scalar*    Value;
    
    bool operator==(const SystemEntry& e) {
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
    
    
    //map functions. Note that you can't map to a jacobi. This is for the case that we 
    //switch to sparse matrices later.                
               
    std::vector<SystemEntry<Kernel>> mapParameter(Scalar*& s) {
        s = &m_parameters(++m_parameterOffset);
        std::vector<SystemEntry<Kernel>> result(1);
        result[0] = {m_parameterOffset, s};
        return result;
    };

    template<typename Derived>
    std::vector<SystemEntry<Kernel>> mapParameter(Eigen::Map<Derived>& map) {
        new(&map) Eigen::Map<Derived>(&m_parameters(++m_parameterOffset));
        
        std::vector<SystemEntry<Kernel>> result(map.rows());
        for (int i = 0; i < map.rows(); ++i)
            result[i] = {m_parameterOffset + i, &m_parameters(m_parameterOffset + i)};

        m_parameterOffset+= 2;
                
        return result;
    };
    
    SystemEntry<Kernel> mapResidual(Scalar* s) {
        s = &m_residuals(++m_residualOffset);
        return {m_parameterOffset, s};
    };
    
    void registerJacobiValue(int row, int col, Scalar s) {
        s = &m_jacobi(row, col);
    };    
    
    void setupJacobi() {
        m_jacobi.setZero();
    };
    
    Scalar& jacobiAt(int row, int col) {
       return m_jacobi(row, col);  
    };
    
    //access the vectors and matrices
    VectorX parameter() {return m_parameters;};
    VectorX residuals() {return m_residuals;};
    MatrixX jacobi()    {return m_jacobi;};    
    
private:
    int m_parameterCount, m_equationCount;
    int m_parameterOffset = -1, m_residualOffset  = -1;
    VectorX m_parameters;
    VectorX m_residuals;
    MatrixX m_jacobi;
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






