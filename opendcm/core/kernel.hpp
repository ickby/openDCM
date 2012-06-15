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

#ifndef DCM_KERNEL_H
#define DCM_KERNEL_H

#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Geometry>

#include <iostream>


namespace dcm {

namespace E = Eigen;


template<typename Kernel>
struct Dogleg {

    typedef typename Kernel::number_type number_type;

    bool solve(typename Kernel::MappedEquationSystem& sys) {

        if(!sys.isValid()) return false;

        number_type tolg=1e-80, tolx=1e-80, tolf=1e-10;

        typename Kernel::Vector g(sys.m_params), h_sd(sys.m_params),
                 h_gn(sys.m_params), h_dl(sys.m_params), F_old(sys.m_eqns);
        typename Kernel::Matrix J_old(sys.m_eqns, sys.m_params);

        sys.recalculate();
	
	std::stringstream stream;
	      stream<<"parameter: "<<std::endl<<sys.Parameter.transpose()<<std::endl<<std::endl;
 	      Base::Console().Message("%s", stream.str().c_str());

        number_type err = sys.Residual.norm();

        F_old = sys.Residual;
        J_old = sys.Jacobi;
        g = sys.Jacobi.transpose()*(sys.Residual);

        // get the infinity norm fx_inf and g_inf
        number_type g_inf = g.template lpNorm<E::Infinity>();
        number_type fx_inf = sys.Residual.template lpNorm<E::Infinity>();

        int maxIterNumber = 1000;//MaxIterations * xsize;
        number_type diverging_lim = 1e6*err + 1e12;

        number_type delta=1;
        number_type alpha=0.;
        number_type nu=2.;
        int iter=0, stop=0, reduce=0;
        while(!stop) {

            // check if finished
            if(fx_inf <= tolf)  // Success
                stop = 1;
            else if(g_inf <= tolg)
                stop = 2;
            else if(delta <= tolx*(tolx + sys.Parameter.norm()))
                stop = 3;
            else if(iter >= maxIterNumber)
                stop = 4;
            else if(err > diverging_lim || err != err) {  // check for diverging and NaN
                stop = 6;
            } else {

                // get the steepest descent stepsize and direction
                alpha = g.squaredNorm()/(sys.Jacobi*g).squaredNorm();
                h_sd  = -g;

                // get the gauss-newton step
                h_gn = sys.Jacobi.fullPivLu().solve(-sys.Residual);

                // compute the dogleg step
                if(h_gn.norm() <= delta) {
		  // std::cout<<"Gauss Newton"<<std::endl;
                    h_dl = h_gn;
                    if(h_dl.norm() <= tolx*(tolx + sys.Parameter.norm())) {
                        stop = 5;
                        break;
                    }
                } else if((alpha*h_sd).norm() >= delta) {
		  // std::cout<<"Steepest descent"<<std::endl;
                    //h_dl = alpha*h_sd;
                    h_dl = (delta/(h_sd.norm()))*h_sd;
                    //die theorie zu dogleg sagt: h_dl = (delta/(h_sd.norm()))*h_sd;
                    //wir gehen aber den klassichen steepest descent weg
                } else {
		  // std::cout<<"Dogleg"<<std::endl;
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
            }

            // see if we are already finished
            if(stop)
                break;

            // calculate the linear model
            number_type dL = 0.5*std::pow(err,2) - 0.5*(sys.Residual + sys.Jacobi*h_dl).squaredNorm();

            // get the new values
            sys.Parameter += h_dl;
            sys.recalculate();
	    
	    // std::cout<<"Parameter Transposed:"<<std::endl<<sys.Parameter.transpose()<<std::endl;
	    // std::cout<<"Residual:"<<std::endl<<sys.Residual<<std::endl;
	    // std::cout<<"Jacobi:"<<std::endl<<sys.Jacobi<<std::endl;
	    
            //calculate the update ratio
            number_type err_new = sys.Residual.norm();
            number_type dF = err - err_new;
            number_type rho = dF/dL;
	    //std::cout<<"rho: "<<rho<<std::endl;

            if(dF > 0 && dL > 0) {

                F_old = sys.Residual;
                J_old = sys.Jacobi;
		
	      std::stringstream stream;
	      stream<<"jacobi: "<<std::endl<<J_old<<std::endl<<"residual:"<<std::endl<<F_old<<std::endl;
	      stream<<"update: "<<std::endl<<h_dl.transpose()<<std::endl;
	      stream<<"delta: "<<delta<<std::endl<<std::endl;
	      
	      Base::Console().Message("%s", stream.str().c_str());


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
                rho = -1;
            }

            // update delta
            if(rho>0.75) {
                delta = std::max(delta,3*h_dl.norm());
                //delta = 3*delta;
                nu = 2;
            } else if(rho < 0.25) {
                delta = delta/nu;
                nu = 2*nu;
            }
            // std::cout<<"Delta: "<<delta<<std::endl<<std::endl;
            // count this iteration and start again
            iter++;
        }
        // std::cout<<"Iterations used: "<<iter<<std::endl<<std::endl;
        Base::Console().Message("residual: %e, reason: %d, iterations: %d\n", err, stop, iter);
        
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
        int m_params, m_eqns;

        MappedEquationSystem(int p, int e) : Jacobi(e, p),
            Parameter(p), Residual(e), m_params(p), m_eqns(e) {};

        void setParameterMap(int offset, int number, VectorMap& map) {
            new(&map) VectorMap(&Parameter(offset), number, DynStride(1,1));
        };
        void setParameterMap(int offset, Vector3Map& map) {
            new(&map) Vector3Map(&Parameter(offset));
        };
        void setResidualMap(int eqn, VectorMap& map) {
            new(&map) VectorMap(&Residual(eqn), 1, DynStride(1,1));
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
        return ((t1-t2) < 0.001);
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

#endif //DCM_KERNEL_H
