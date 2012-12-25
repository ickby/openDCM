/*
    openDCM, dimensional constraint manager
    Copyright (C) 2012  Stefan Troeger <stefantroeger@gmx.net>

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/

#ifndef GCM_KERNEL_H
#define GCM_KERNEL_H

#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Geometry>

#include <iostream>
//#include <../FreeCAD/src/Base/Console.h>

#include <boost/math/special_functions/fpclassify.hpp>
#include <time.h>

namespace dcm {

namespace E = Eigen;

template<typename Kernel>
struct Dogleg {

    typedef typename Kernel::number_type number_type;
    number_type tolg, tolx, tolf;

    Dogleg() : tolg(1e-80), tolx(1e-10), tolf(1e-5) {};

    template <typename Derived, typename Derived2, typename Derived3, typename Derived4>
    int calculateStep(const Eigen::MatrixBase<Derived>& g, const Eigen::MatrixBase<Derived3>& jacobi,
                      const Eigen::MatrixBase<Derived4>& residual, Eigen::MatrixBase<Derived2>& h_dl,
                      const double delta) {

        // get the steepest descent stepsize and direction
        const double alpha(g.squaredNorm()/(jacobi*g).squaredNorm());
        const typename Kernel::Vector h_sd  = -g;

        // get the gauss-newton step
        const typename Kernel::Vector h_gn = (jacobi).fullPivLu().solve(-residual);

        // compute the dogleg step
        if(h_gn.norm() <= delta) {
            h_dl = h_gn;
        } else if((alpha*h_sd).norm() >= delta) {
            //h_dl = alpha*h_sd;
            h_dl = (delta/(h_sd.norm()))*h_sd;
        } else {
            //compute beta
            number_type beta = 0;
            typename Kernel::Vector a = alpha*h_sd;
            typename Kernel::Vector b = h_gn;
            number_type c = a.transpose()*(b-a);
            number_type bas = (b-a).squaredNorm(), as = a.squaredNorm();
            if(c<0) {
                beta = -c+std::sqrt(std::pow(c,2)+bas*(std::pow(delta,2)-as));
                beta /= bas;
            } else {
                beta = std::pow(delta,2)-as;
                beta /= c+std::sqrt(std::pow(c,2) + bas*(std::pow(delta,2)-as));
            };

            // and update h_dl and dL with beta
            h_dl = alpha*h_sd + beta*(b-a);
        }
        return 0;
    };

    bool solve(typename Kernel::MappedEquationSystem& sys) {
        //std::cout<<"start solving"<<std::endl;
        clock_t start = clock();
        clock_t inc_rec = clock();

        if(!sys.isValid()) return false;

        bool translate = true;

        //Base::Console().Message("\nparams: %d, rot_params: %d, trans_params: %d\n", sys.m_params, sys.m_rot_params, sys.m_trans_params);
        typename Kernel::Vector h_dl, F_old(sys.m_eqns), g(sys.m_eqns);
        typename Kernel::Matrix J_old(sys.m_eqns, sys.m_params);

        sys.recalculate();

        std::stringstream s;
        s<<"initial jacobi: "<<std::endl<<sys.Jacobi<<std::endl;
        s<<"residual: "<<sys.Residual.transpose()<<std::endl;
        // Base::Console().Message("%s", s.str().c_str());

        number_type err = sys.Residual.norm();

        F_old = sys.Residual;
        J_old = sys.Jacobi;

        g = sys.Jacobi.transpose()*(sys.Residual);

        // get the infinity norm fx_inf and g_inf
        number_type g_inf = g.template lpNorm<E::Infinity>();
        number_type fx_inf = sys.Residual.template lpNorm<E::Infinity>();

        int maxIterNumber = 10000;//MaxIterations * xsize;
        number_type diverging_lim = 1e6*err + 1e12;

        number_type delta=5;
        number_type nu=2.;
        int iter=0, stop=0, reduce=0, unused=0;


        while(!stop) {

            // check if finished
            if(fx_inf <= tolf*sys.Scaling)  // Success
                stop = 1;
            else if(g_inf <= tolg)
                stop = 2;
            else if(delta <= tolx)
                stop = 3;
            else if(iter >= maxIterNumber)
                stop = 4;
            else if(err > diverging_lim || err != err) {  // check for diverging and NaN
                stop = 6;
            }

            // see if we are already finished
            if(stop)
                break;

            number_type err_new;
            number_type dF=0, dL=0;
            number_type rho;

            //get the update step
            calculateStep(g, sys.Jacobi,  sys.Residual, h_dl, delta);

            // calculate the linear model
            dL = 0.5*sys.Residual.norm() - 0.5*(sys.Residual + sys.Jacobi*h_dl).norm();

            // get the new values
            sys.Parameter += h_dl;

            clock_t start_rec = clock();
            sys.recalculate();
            clock_t end_rec = clock();
            inc_rec += end_rec-start_rec;

            //calculate the translation update ratio
            err_new = sys.Residual.norm();
            dF = err - err_new;
            rho = dF/dL;

            if(dF<=0 || dL<=0)  rho = -1;
            // update delta
            if(rho>0.75) {
                delta = std::max(delta,3*h_dl.norm());
                nu = 2;
            } else if(rho < 0.25) {
                delta = delta/nu;
                nu = 2*nu;
            }

            //Base::Console().Message("delta: %e, error: %e\n", delta, err);


            if(dF > 0 && dL > 0) {

                F_old = sys.Residual;
                J_old = sys.Jacobi;

                err = err_new;

                g = sys.Jacobi.transpose()*(sys.Residual);

                // get infinity norms
                g_inf = g.template lpNorm<E::Infinity>();
                fx_inf = sys.Residual.template lpNorm<E::Infinity>();

                //stream<<"accepted, dr and dt:"<<delta_r<<", "<<delta_t<<std::endl<<std::endl;

//                 std::stringstream s;
//                 s<<"accepted jacobi: "<<std::endl<<sys.Jacobi<<std::endl;
// 		s<<"step: "<<h_dl.transpose()<<std::endl;
//                 s<<"residual: "<<sys.Residual.transpose()<<std::endl;
// 		s<<"delta: "<<delta;
// 		Base::Console().Message("%s", s.str().c_str());
                // count this iteration and start again

            } else {
                // std::cout<<"Step Rejected"<<std::endl;
                sys.Residual = F_old;
                sys.Jacobi = J_old;
                sys.Parameter -= h_dl;
                unused++;
            }

            iter++;
        }

        std::stringstream ss;
        ss<<"end jacobi: "<<std::endl<<sys.Jacobi<<std::endl;
        //Base::Console().Message("%s", ss.str().c_str());

        clock_t end = clock();
        double ms = (double(end-start) * 1000.) / double(CLOCKS_PER_SEC);
        double ms_rec = (double(inc_rec-start) * 1000.) / double(CLOCKS_PER_SEC);
        //Base::Console().Message("residual: %e, reason: %d, iterations: %d, time in ms: %f, recalc time %f\n",
        //                       err, stop, iter, ms, ms_rec);
        //std::cout<<"DONE solving"<<std::endl;

        if(stop == 1) return true;
        return false; //TODO:throw
    }
};

template<typename Scalar, template<class> class Solver = Dogleg>
struct Kernel {

    //basics
    typedef Scalar number_type;

    //Linear algebra types
    typedef E::Matrix<Scalar, 3, 1> Vector3;
    typedef E::Matrix<Scalar, 1, 3> CVector3;
    typedef E::Matrix<Scalar, 3, 3> Matrix3;
    typedef E::Matrix<Scalar, E::Dynamic, 1> Vector;
    typedef E::Matrix<Scalar, 1, E::Dynamic> CVector;
    typedef E::Matrix<Scalar, E::Dynamic, E::Dynamic> Matrix;

    //mapped types
    typedef E::Stride<E::Dynamic, E::Dynamic> DynStride;
    typedef E::Map< Vector3 > Vector3Map;
    typedef E::Map< CVector3> CVector3Map;
    typedef E::Map< Matrix3 > Matrix3Map;
    typedef E::Map< Vector, 0, DynStride > VectorMap;
    typedef E::Map< CVector, 0, DynStride > CVectorMap;
    typedef E::Map< Matrix, 0, DynStride > MatrixMap;

    //Special types
    typedef E::Quaternion<Scalar>   Quaternion;
    typedef E::Matrix<Scalar, 3, 9> Matrix39;
    typedef E::Map< Matrix39 >      Matrix39Map;
    typedef E::Block<Matrix>	    MatrixBlock;

    struct MappedEquationSystem {

        Matrix Jacobi;
        Vector Parameter;
        Vector Residual;
        number_type Scaling;
        int m_params, m_eqns; //total amount
        int m_param_offset, m_eqn_offset;   //current positions while creation

        MappedEquationSystem(int params, int equations)
            : Jacobi(equations, params),
              Parameter(params), Residual(equations),
              m_params(params), m_eqns(equations) {

            m_param_offset = 0;
            m_eqn_offset = 0;

            Jacobi.setZero(); //important as some places are never written
        };

        int setParameterMap(int number, VectorMap& map) {

            new(&map) VectorMap(&Parameter(m_param_offset), number, DynStride(1,1));
            m_param_offset += number;
            return m_param_offset-number;
        };
        int setParameterMap(Vector3Map& map) {

            new(&map) Vector3Map(&Parameter(m_param_offset));
            m_param_offset += 3;
            return m_param_offset-3;
        };
        int setResidualMap(VectorMap& map) {
            new(&map) VectorMap(&Residual(m_eqn_offset), 1, DynStride(1,1));
            return m_eqn_offset++;
        };
        void setJacobiMap(int eqn, int offset, int number, CVectorMap& map) {
            new(&map) CVectorMap(&Jacobi(eqn, offset), number, DynStride(0,m_eqns));
        };
        void setJacobiMap(int eqn, int offset, int number, VectorMap& map) {
            new(&map) VectorMap(&Jacobi(eqn, offset), number, DynStride(0,m_eqns));
        };

        bool isValid() {
            if(!m_params || !m_eqns) return false;
            return true;
        };

        virtual void recalculate() = 0;

    };

    Kernel()  {};

    template <typename DerivedA,typename DerivedB>
    static bool isSame(const E::MatrixBase<DerivedA>& p1,const E::MatrixBase<DerivedB>& p2) {
        return ((p1-p2).squaredNorm() < 0.001);
    }
    static bool isSame(number_type t1, number_type t2) {
        return (std::abs(t1-t2) < 0.001);
    }
    template <typename DerivedA,typename DerivedB>
    static bool isOpposite(const E::MatrixBase<DerivedA>& p1,const E::MatrixBase<DerivedB>& p2) {
        return ((p1+p2).squaredNorm() < 0.001);
    }

    static bool solve(MappedEquationSystem& mes) {
        return Solver< Kernel<Scalar, Solver> >().solve(mes);
    };

};

}

#endif //GCM_KERNEL_H


