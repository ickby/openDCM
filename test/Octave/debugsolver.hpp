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

#ifndef GCM_DEBUG_SOLVER_H
#define GCM_DEBUG_SOLVER_H

#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Geometry>

#include <iostream>
#include <fstream>


namespace dcm {

template<typename Kernel>
struct DebugSolver {

    typedef typename Kernel::number_type number_type;
    number_type tolg, tolx, tolf;
#ifdef USE_LOGGING
    src::logger log;
#endif

    DebugSolver() : tolg(1e-80), tolx(1e-10), tolf(1e-5) {
#ifdef USE_LOGGING
        log.add_attribute("Tag", attrs::constant< std::string >("Dogleg"));
#endif
    };

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
        typename Kernel::Vector h_dl, F_old(sys.m_eqns), g(sys.m_eqns), Original(sys.m_params);
        typename Kernel::Matrix J_old(sys.m_eqns, sys.m_params);

        sys.recalculate();

        Original = sys.Parameter;

#ifdef USE_LOGGING
        BOOST_LOG(log)<< "initial jacobi: "<<std::endl<<sys.Jacobi<<std::endl
                      << "residual: "<<sys.Residual.transpose();
#endif

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

            if(dF > 0 && dL > 0) {

                F_old = sys.Residual;
                J_old = sys.Jacobi;

                err = err_new;

                g = sys.Jacobi.transpose()*(sys.Residual);

                // get infinity norms
                g_inf = g.template lpNorm<E::Infinity>();
                fx_inf = sys.Residual.template lpNorm<E::Infinity>();

            } else {
                // std::cout<<"Step Rejected"<<std::endl;
                sys.Residual = F_old;
                sys.Jacobi = J_old;
                sys.Parameter -= h_dl;
                unused++;
            }

            iter++;
        }
        clock_t end = clock();
        double ms = (double(end-start) * 1000.) / double(CLOCKS_PER_SEC);

#ifdef USE_LOGGING
        BOOST_LOG(log) <<"Done solving: "<<err<<", iter: "<<iter;
#endif

        if(/*stop != */1) {

            int jcount = 1000;

            std::ofstream str("/home/stefan/Projects/openDCM/test/Octave/output.m");
            if(!str.is_open())
                std::cout<<"file not opend!"<<std::endl;

            //Reset and start varying all parameter
            sys.Parameter = Original;
            for(int i=0; i<Original.rows(); i++) {

                typename Kernel::Vector param_gradient(jcount);
                typename Kernel::Vector result(jcount);

                str<<"parameter("<<i+1<<",:) = [";
                //varying parameter i
                for(int j=0; j<jcount; j++) {
                    sys.Parameter(i) = Original(i) + double(jcount/2-j)* 1./double(jcount);
                    str<<sys.Parameter(i)<<" ";
                    //calculate the gradient
                    sys.recalculate();
                    result(j) = sys.Residual(0);
                    param_gradient(j) = sys.Jacobi(0,i);
                };
                str<<"];"<<std::endl;
                //write the results
                str<<"results("<<i+1<<",:) = [";
                for(int j=0; j<jcount; j++)
                    str<<result(j)<<" ";
                str<<"];"<<std::endl;
                //write the gradients
                str<<"gradient("<<i+1<<",:) = [";
                for(int j=0; j<jcount; j++)
                    str<<param_gradient(j)<<" ";
                str<<"];"<<std::endl;

                str<<std::endl;

            };

            str.close();
        };

        if(stop == 1) return true;
        return false; //TODO:throw
    }
};


}

#endif //GCM_KERNEL_H
