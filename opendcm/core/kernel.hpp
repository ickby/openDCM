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

#include "transformation.hpp"
#include "logging.hpp"

namespace dcm {
namespace numeric {
    
template< typename Kernel >
struct Parameter {
    
    int                         Index;
    typename Kernel::Scalar*    Value;
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
struct Eigen3Kernel {

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






