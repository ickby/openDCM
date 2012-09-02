/*
    openGCM, geometric constraint manager
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


namespace gcm {

namespace E = Eigen;

enum ParameterType {
    Rotation,
    Translation,
    Anything
};

template<typename Kernel>
struct Dogleg {

    typedef typename Kernel::number_type number_type;
    number_type tolg, tolx, tolf;

    Dogleg() : tolg(1e-80), tolx(1e-10), tolf(1e-10) {};

    template <typename Derived, typename Derived2, typename Derived3, typename Derived4>
    int calculateStep(const Eigen::MatrixBase<Derived>& g, const Eigen::MatrixBase<Derived3>& jacobi,
                      const Eigen::MatrixBase<Derived4>& residual, Eigen::MatrixBase<Derived2>& h_dl,
                      const double delta) {

        // std::stringstream stream;
        // stream<</*"g: "<<std::endl<<g<<std::endl<<*/"jacobi: "<<std::endl<<jacobi<<std::endl;
        // stream<<"residual: "<<std::endl<<residual<<std::endl;


        // get the steepest descent stepsize and direction
        const double alpha(g.squaredNorm()/(jacobi*g).squaredNorm());
        const typename Kernel::Vector h_sd  = -g;

        // get the gauss-newton step
        const typename Kernel::Vector h_gn = (jacobi).fullPivLu().solve(-residual);

        // compute the dogleg step
        if(h_gn.norm() <= delta) {
            // std::cout<<"Gauss Newton"<<std::endl;
            h_dl = h_gn;
            //if(h_dl.norm() <= tolx*(tolx + sys.Parameter.norm())) {
            //    return 5;
            //}
        } else if((alpha*h_sd).norm() >= delta) {
            // std::cout<<"Steepest descent"<<std::endl;
            //h_dl = alpha*h_sd;
            h_dl = (delta/(h_sd.norm()))*h_sd;
            //stream<<"h_dl steepest: "<<std::endl<<h_dl<<std::endl;
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
        // stream<<"jacobi*h_dl"<<std::endl<<(jacobi*h_dl)<<std::endl<<std::endl;
        //stream<<"h_dl:"<<std::endl<<h_dl<<std::endl<<std::endl;
        //  Base::Console().Message("%s", stream.str().c_str());
        return 0;
    };

    bool solve(typename Kernel::MappedEquationSystem& sys) {
        //std::cout<<"start solving"<<std::endl;
        clock_t start = clock();


        if(!sys.isValid()) return false;

        int npt = sys.m_params+sys.m_trans_params;
        int npr = sys.m_rot_params;
        bool translate = true;

        //Base::Console().Message("\nparams: %d, rot_params: %d, trans_params: %d\n", sys.m_params, sys.m_rot_params, sys.m_trans_params);
        typename Kernel::Vector h_dlt(npt), h_dlr(npr),
                 F_old(sys.m_eqns), g(sys.m_eqns), BFR_Res(sys.m_eqns);
        typename Kernel::Matrix J_old(sys.m_eqns, npt+npr), BFR_J(sys.m_eqns, npt+npr);

        sys.recalculate();

        // std::stringstream stream;
        // stream<<"start jacobi: "<<std::endl<<sys.Jacobi<<std::endl;
        // stream<<"parameter: "<<std::endl<<sys.Parameter.transpose()<<std::endl<<std::endl;
        //Base::Console().Message("%s", stream.str().c_str());
        // stream.str(std::string());

        number_type err = sys.Residual.norm();

        F_old = sys.Residual;
        J_old = sys.Jacobi;

        g = sys.Jacobi.transpose()*(sys.Residual);
        if((npt != 0) && (g.head(npt).norm() == 0)) translate=false;   //means that translations are existend but have no influence

        // get the infinity norm fx_inf and g_inf
        number_type g_inf = g.template lpNorm<E::Infinity>();
        number_type fx_inf = sys.Residual.template lpNorm<E::Infinity>();

        int maxIterNumber = 1000;//MaxIterations * xsize;
        number_type diverging_lim = 1e6*err + 1e12;

        number_type delta_t=5, delta_r=0.1;
        number_type nu_t=2., nu_r=2.;
        int iter=0, stop=0, reduce=0, unused_lin=0, unused_rot=0;

//                 std::stringstream stream;
//                 stream<<"init jacobi:"<<std::endl<<sys.Jacobi<<std::endl<<std::endl;
//         	stream<<"init parameter:"<<std::endl<<sys.Parameter.transpose()<<std::endl<<std::endl;
//         //        stream<<std::fixed<<std::setprecision(5)<<"delta_t: "<<delta_t<<",   delta_r: " << delta_r;
//                 stream<<"init residual:"<<sys.Residual.transpose()<<std::endl;
//                 Base::Console().Message("%s", stream.str().c_str());

        while(!stop) {

            // check if finished
            if(fx_inf <= tolf)  // Success
                stop = 1;
            else if(g_inf <= tolg)
                stop = 2;
            else if(delta_t <= tolx && delta_r <= tolx)
                stop = 3;
            else if(iter >= maxIterNumber)
                stop = 4;
            else if(err > diverging_lim || err != err) {  // check for diverging and NaN
                stop = 6;
            }

            // see if we are already finished
            if(stop)
                break;

            std::stringstream stream, stream2;

            number_type err_new_r;
            if(npr) {
                //get the update step
                calculateStep(g.tail(npr), sys.Jacobi.block(0, npt, sys.m_eqns, npr), sys.Residual, h_dlr, delta_r);

                //calculate linear model
                BFR_Res = sys.Residual;
                BFR_J   = sys.Jacobi;


                //stream<<"Jacobi npt: "<<std::endl<<sys.Jacobi.block(0, npt, sys.m_eqns, npr)<<std::endl;
                //stream<<"Update: "<<std::endl<<h_dlr<<std::endl;

                sys.Parameter.tail(npr) += h_dlr;
                sys.recalculate();

                //calculate the rotation update ratio
                err_new_r = sys.Residual.norm();
            } else {
                err_new_r = err;
                h_dlr.setZero();
            }

            //stream<<"after npr residual: "<<sys.Residual.transpose()<<std::endl;


            number_type dF_t=0, dL_t=0;
            number_type rho_t, err_new_t;
            if(npt && translate) {
                //get the update step
                calculateStep(g.head(npt), sys.Jacobi.block(0, 0, sys.m_eqns, npt),
                              sys.Residual, h_dlt, delta_t);

                // calculate the linear model
                dL_t = 0.5*sys.Residual.norm() - 0.5*(sys.Residual + sys.Jacobi.block(0, 0, sys.m_eqns, npt)*h_dlt).norm();

                std::stringstream stream;
                //stream<<"residual:"<<std::endl<<sys.Residual<<std::endl<<std::endl;
                //stream<<"residual + jacobi*h_dl"<<std::endl<<(sys.Residual + sys.Jacobi.block(0, 0, sys.m_eqns, npt)*h_dlt)<<std::endl<<std::endl;
                //Base::Console().Message("%s", stream.str().c_str());

                // get the new values
                sys.Parameter.head(npt) += h_dlt;
                sys.recalculate();

                //calculate the translation update ratio
                err_new_t = sys.Residual.norm();
                dF_t = err_new_r - err_new_t;
                rho_t = dF_t/dL_t;

                if(dF_t<=0 || dL_t<=0) {
                    rho_t = -1;
                    unused_lin++;
                }
                // update delta
                if(rho_t>0.75) {
                    delta_t = std::max(delta_t,3*h_dlt.norm());
                    nu_t = 2;
                } else if(rho_t < 0.25) {
                    delta_t = delta_t/nu_t;
                    nu_t = 2*nu_t;
                }
            } else {
                err_new_t = err_new_r;
                h_dlt.setZero();
            }

            //stream<<"after npt residual: "<<sys.Residual.transpose()<<std::endl;

            number_type dL, dF;
            typename Kernel::Vector h_dl(npt+npr);
            h_dl.head(npt) = h_dlt;
            h_dl.tail(npr) = h_dlr;
            if(npr) {
                // calculate the linear model
                dL = 0.5*BFR_Res.norm() - 0.5*(BFR_Res + BFR_J*h_dl).norm();


                //calculate the update ratio
                dF = err - err_new_t;
                number_type rho = dF/dL;

                if(dF<=0 || dL<=0) {
                    rho = -1;
                    unused_rot++;
                }
                // update delta
                if(rho>0.75) {
                    delta_r = std::max(delta_r,2*h_dlr.norm());
                    nu_r = 2;
                } else if(rho < 0.25) {
                    delta_r = delta_r/nu_r;
                    nu_r = 2*nu_r;
                }
            } else {
                dF = dF_t;
                dL = dL_t;
            }

            if(dF > 0 && dL > 0) {

                F_old = sys.Residual;
                J_old = sys.Jacobi;

                err = err_new_t;

                g = sys.Jacobi.transpose()*(sys.Residual);

                // get infinity norms
                g_inf = g.template lpNorm<E::Infinity>();
                fx_inf = sys.Residual.template lpNorm<E::Infinity>();

                //stream<<"accepted, dr and dt:"<<delta_r<<", "<<delta_t<<std::endl<<std::endl;

                std::stringstream s;
                s<<"accepted jacobi: "<<std::endl<<sys.Jacobi<<std::endl;
                s<<"residual: "<<sys.Residual.transpose()<<std::endl;
                Base::Console().Message("%s", s.str().c_str());
                // count this iteration and start again

            } else {
                // std::cout<<"Step Rejected"<<std::endl;
                sys.Residual = F_old;
                sys.Jacobi = J_old;
                sys.Parameter -= h_dl;
            }


            // std::stringstream stream;
            // stream<<std::fixed<<std::setprecision(5)<<"delta_t: "<<delta_t<<",   delta_r: " << delta_r;
            // stream<<",  parameter:"<<sys.Parameter.transpose()<<std::endl;
            //Base::Console().Message("%s", stream.str().c_str());
            // std::cout<<"Delta: "<<delta<<std::endl<<std::endl;
            // count this iteration and start again
            iter++;
        }
//         std::stringstream s;
//         s<<"end jacobi: "<<std::endl<<sys.Jacobi<<std::endl;
//         s<<"parameter: "<<std::endl<<sys.Parameter.transpose()<<std::endl;
//         s<<"residual: "<<std::endl<<sys.Residual<<std::endl<<std::endl;
//         Base::Console().Message("%s", s.str().c_str());
        // std::cout<<"Iterations used: "<<iter<<std::endl<<std::endl;
        clock_t end = clock();
        double ms = (double(end-start) * 1000.) / double(CLOCKS_PER_SEC);
        Base::Console().Message("residual: %e, reason: %d, iterations: %d, time in ms: %f, unused lin: %d, unused rot: %d\n",
                                err, stop, iter, ms, unused_lin, unused_rot);
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
	Vector Scaling;
        int m_params, m_rot_params, m_trans_params, m_eqns; //total amount
        int m_rot_offset, m_trans_offset, m_param_offset, m_eqn_offset;   //current positions while creation
        number_type m_scale;

        MappedEquationSystem(int params, int rotparams, int transparams, int equations)
            : Jacobi(equations, params+rotparams+transparams), Scaling(equations),
              Parameter(params+rotparams+transparams), Residual(equations),
              m_params(params), m_rot_params(rotparams), m_trans_params(transparams), m_eqns(equations) {
            m_rot_offset = params+transparams;
            m_trans_offset = params;
            m_param_offset = 0;
            m_eqn_offset = 0;

            Jacobi.setZero(); //important as some places are never written
        };

        int setParameterMap(ParameterType t, int number, VectorMap& map) {
            if(t==Anything) {
                new(&map) VectorMap(&Parameter(m_param_offset), number, DynStride(1,1));
                m_param_offset += number;
                return m_param_offset-number;
            } else if(t==Rotation) {
                new(&map) VectorMap(&Parameter(m_rot_offset), number, DynStride(1,1));
                m_rot_offset += number;
                return m_rot_offset-number;
            } else if(t==Translation) {
                new(&map) VectorMap(&Parameter(m_trans_offset), number, DynStride(1,1));
                m_trans_offset += number;
                return m_trans_offset-number;
            }
        };
        int setParameterMap(ParameterType t, Vector3Map& map) {
            if(t==Anything) {
                new(&map) Vector3Map(&Parameter(m_param_offset));
                m_param_offset += 3;
                return m_param_offset-3;
            } else if(t==Rotation) {
                new(&map) Vector3Map(&Parameter(m_rot_offset));
                m_rot_offset += 3;
                return m_rot_offset-3;
            } else if(t==Translation) {
                new(&map) Vector3Map(&Parameter(m_trans_offset));
                m_trans_offset += 3;
                return m_trans_offset-3;
            }
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
            if(!(m_params+m_trans_params+m_rot_params) || !m_eqns) return false;

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


template<typename Kernel>
struct Simple {

    typedef typename Kernel::number_type number_type;
    number_type tolg, tolx, tolf;

    Simple() : tolg(1e-80), tolx(1e-10), tolf(1e-10) {};

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


        if(!sys.isValid()) return false;

        int npt = sys.m_params+sys.m_trans_params;
        int npr = sys.m_rot_params;
        bool translate = true;

        //Base::Console().Message("\nparams: %d, rot_params: %d, trans_params: %d\n", sys.m_params, sys.m_rot_params, sys.m_trans_params);
        typename Kernel::Vector h_dl, F_old(sys.m_eqns), g(sys.m_eqns);
        typename Kernel::Matrix J_old(sys.m_eqns, npt+npr);

        sys.recalculate();
	
 	std::stringstream s;
                 s<<"initial jacobi: "<<std::endl<<sys.Jacobi<<std::endl;
                 s<<"residual: "<<sys.Residual.transpose()<<std::endl;
                 Base::Console().Message("%s", s.str().c_str());

        number_type err = sys.Residual.norm();

        F_old = sys.Residual;
        J_old = sys.Jacobi;

        g = sys.Jacobi.transpose()*(sys.Residual);

        // get the infinity norm fx_inf and g_inf
        number_type g_inf = g.template lpNorm<E::Infinity>();
        number_type fx_inf = sys.Residual.template lpNorm<E::Infinity>();

        int maxIterNumber = 1000;//MaxIterations * xsize;
        number_type diverging_lim = 1e6*err + 1e12;

        number_type delta=5;
        number_type nu=2.;
        int iter=0, stop=0, reduce=0, unused=0;


        while(!stop) {

            // check if finished
            if(fx_inf <= tolf)  // Success
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
            sys.recalculate();

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
            
            Base::Console().Message("delta: %e, error: %e\n", delta, err);


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
	clock_t end = clock();
        double ms = (double(end-start) * 1000.) / double(CLOCKS_PER_SEC);
        Base::Console().Message("residual: %e, reason: %d, iterations: %d, time in ms: %f, unused %d\n",
                                err, stop, iter, ms, unused);
        //std::cout<<"DONE solving"<<std::endl;

        if(stop == 1) return true;
        return false; //TODO:throw
    }
};

template<typename Kernel>
struct SimpleScaled {

    typedef typename Kernel::number_type number_type;
    number_type tolg, tolx, tolf;

    SimpleScaled() : tolg(1e-80), tolx(1e-10), tolf(1e-10) {};

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


        if(!sys.isValid()) return false;

        int npt = sys.m_params+sys.m_trans_params;
        int npr = sys.m_rot_params;
        bool translate = true;

        //Base::Console().Message("\nparams: %d, rot_params: %d, trans_params: %d\n", sys.m_params, sys.m_rot_params, sys.m_trans_params);
        typename Kernel::Vector h_dl, F_old(sys.m_eqns), g(sys.m_eqns);
        typename Kernel::Matrix J_old(sys.m_eqns, npt+npr);
	
	typename Kernel::Matrix rmatscale(sys.m_eqns, npr+npt);
        for(int i=0; i<(npr+npt); i++)
            rmatscale.col(i) = sys.Scaling;
	

        sys.recalculate();
	sys.Jacobi   = sys.Jacobi.cwiseProduct(rmatscale);
	sys.Residual = sys.Residual.cwiseProduct(sys.Scaling);
	
	
 	std::stringstream s;
        s<<"initial jacobi: "<<std::endl<<sys.Jacobi<<std::endl;
        s<<"residual: "<<sys.Residual.transpose()<<std::endl;
        Base::Console().Message("%s", s.str().c_str());

        number_type err = sys.Residual.norm();

        F_old = sys.Residual;
        J_old = sys.Jacobi;

        g = sys.Jacobi.transpose()*(sys.Residual);

        // get the infinity norm fx_inf and g_inf
        number_type g_inf = g.template lpNorm<E::Infinity>();
        number_type fx_inf = sys.Residual.template lpNorm<E::Infinity>();

        int maxIterNumber = 1000;//MaxIterations * xsize;
        number_type diverging_lim = 1e6*err + 1e12;

        number_type delta=5;
        number_type nu=2.;
        int iter=0, stop=0, reduce=0, unused=0;


        while(!stop) {

            // check if finished
            if(fx_inf <= tolf)  // Success
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
            sys.recalculate();
	    sys.Jacobi   = sys.Jacobi.cwiseProduct(rmatscale);
	    sys.Residual = sys.Residual.cwiseProduct(sys.Scaling);

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
            
            Base::Console().Message("delta: %e, error: %e\n", delta, err);


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
	clock_t end = clock();
        double ms = (double(end-start) * 1000.) / double(CLOCKS_PER_SEC);
        Base::Console().Message("residual: %e, reason: %d, iterations: %d, time in ms: %f, unused %d\n",
                                err, stop, iter, ms, unused);
	
	                 std::stringstream ss;
//                 s<<"accepted jacobi: "<<std::endl<<sys.Jacobi<<std::endl;
// 		s<<"step: "<<h_dl.transpose()<<std::endl;
                 ss<<"residual: "<<sys.Residual.transpose()<<std::endl;
// 		s<<"delta: "<<delta;
 		Base::Console().Message("%s", ss.str().c_str());
        //std::cout<<"DONE solving"<<std::endl;

        if(stop == 1) return true;
        return false; //TODO:throw
    }
};


}

#endif //GCM_KERNEL_H


