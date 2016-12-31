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
#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/exception/errinfo_errno.hpp>

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
        return (Value == e.Value);
    };
    
    operator Scalar&() {return *Value;};
    
    void operator=(const Scalar& in) {
       *Value = in;
    }
};

template< typename Kernel >
struct MatrixEntry {
    
    typedef typename Kernel::Scalar Scalar;
    
    int     pRow;
    int     Column;
    Scalar* Value;
    
    bool operator==(const MatrixEntry& e) {
        //comparing the value addresses is enough, as if tey point to the same memory they
        //must be the same system entry
        return (Value == Value);
    };
    
    operator Scalar&() {return *Value;};
    
    void operator=(const Scalar& in) {
       *Value = in;
    }
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
    
    //access some stats
    int equationCount()  {return m_equationCount;};
    int parameterCount() {return m_parameterCount;};
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
     * @param sys LinearSystem the calculatable is initalized with
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
    unsigned int m_parameterCount = 0; //how many parameters are added by this equation?

};

/**
 * @brief A sequential vector of calculatables 
 * 
 * This struct is an Executable and works like the shedule vectors. It does however extend their 
 * functionality with convienience methods for \ref Calculatble objects. It allows to initialize 
 * all Calculatable objects in the vector at once.
 * 
 */
template<typename Kernel>
struct CalculatableSequentialVector : public shedule::_SequentialVector<
                                                std::vector<std::shared_ptr<Calculatable<Kernel>>>> {
    
    typedef shedule::_SequentialVector<std::vector<std::shared_ptr<Calculatable<Kernel>>>> Base;
    
    /**
     * @brief Initialization of all Calculatables
     * 
     * Calls the init function of all stored objects     * 
     * @param sys LinearSystem the calculatables are initalized with
     * @return void
     */
    void init(LinearSystem<Kernel>& sys) {
        for(auto calc : Base::m_executables)
            calc->init(sys);
    };
    
    /**
     * @brief Number of free parameters all stored equation need
     * 
     * Return the amount of parameters all storeds equation will need in the solving process.
     * @return void
     */
    unsigned int newParameterCount() {return m_parameterCount;};
    
    void addExecutable(std::shared_ptr<Calculatable<Kernel>> ex) {
        m_parameterCount += ex->newParameterCount();
        Base::m_executables.push_back(ex);
    };
    
private:
    unsigned int m_parameterCount = 0;
};
    
//the standart eigen3 dogleg solver
template<typename Kernel>
struct Dogleg  : shedule::Executable {

#ifdef DCM_USE_LOGGING
    dcm_logger log;
#endif
private:

    typedef typename Kernel::Scalar                               Scalar;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1>              Vector;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> Matrix;
    
    std::shared_ptr<CalculatableSequentialVector<Kernel>> m_equations;
    LinearSystem<Kernel> m_system;
    Scalar tolg, tolx, delta, nu, g_inf, fx_inf, err, time;
    Kernel* m_kernel;
    int iter, stop, reduce, unused, counter;
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> h_dl, F_old, g;
    Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> J_old;

    template <typename Derived, typename Derived2, typename Derived3, typename Derived4>
    void calculateStep(const Eigen::MatrixBase<Derived>& g, const Eigen::MatrixBase<Derived3>& jacobi,
                       const Eigen::MatrixBase<Derived4>& residual, Eigen::MatrixBase<Derived2>& h_dl,
                       const double delta);
        
public:
    Dogleg(std::shared_ptr<CalculatableSequentialVector<Kernel>> vec) : m_system(LinearSystem<Kernel>(vec->newParameterCount(), vec->size())),
                                                        m_equations(vec) {}; 

    virtual void execute();
};

template<typename Kernel>
template <typename Derived, typename Derived2, typename Derived3, typename Derived4>
void Dogleg<Kernel>::calculateStep(const Eigen::MatrixBase<Derived>& g, const Eigen::MatrixBase<Derived3>& jacobi,
                                   const Eigen::MatrixBase<Derived4>& residual, Eigen::MatrixBase<Derived2>& h_dl,
                                   const double delta) {

    // get the steepest descent stepsize and direction
    const double alpha(g.squaredNorm()/(jacobi*g).squaredNorm());
    const Vector h_sd  = -g;

    // get the gauss-newton step
    const Vector h_gn = jacobi.fullPivLu().solve(-residual);
    const double eigen_error = (jacobi*h_gn + residual).norm();
#ifdef USE_LOGGING
    if(eigen_error>1e-9)
        BOOST_LOG_SEV(log, error) << "Linear system solved incorrect: "<<eigen_error;
    
    if(!boost::math::isfinite(h_gn.norm())) {
        BOOST_LOG_SEV(log, error) << "Unnormal gauss-newton detected: "<<h_gn.norm();
    }
    else if(!boost::math::isfinite(h_sd.norm())) {
        BOOST_LOG_SEV(log, error) << "Unnormal steepest descent detected: "<<h_sd.norm();
    }
    else if(!boost::math::isfinite(alpha)) {
        BOOST_LOG_SEV(log, error) << "Unnormal alpha detected: "<<alpha;
    } 
    else {
          BOOST_LOG_SEV(log, iteration) << "Steepest descent: "<<h_sd.transpose();
          BOOST_LOG_SEV(log, iteration) << "Gauss-N  descent: "<<h_gn.transpose();
    };

#endif

    // compute the dogleg step
    if(h_gn.norm() <= delta) {
        h_dl = h_gn;
#ifdef USE_LOGGING
        BOOST_LOG_SEV(log, iteration) << "Gauss-Newton step"; 
#endif
    }
    else if((alpha*h_sd).norm() >= delta) {
        //h_dl = alpha*h_sd;
        h_dl = (delta/(h_sd.norm()))*h_sd;
#ifdef USE_LOGGING
        BOOST_LOG_SEV(log, iteration) << "steepest descent step"; 
        if(!boost::math::isfinite(h_dl.norm())) {
            BOOST_LOG_SEV(log, error) << "Unnormal dogleg descent detected: "<<h_dl.norm();
        }

#endif
    }
    else {
        //compute beta
        typename Kernel::Scalar beta = 0;
        Vector a = alpha*h_sd;
        Vector b = h_gn;
        typename Kernel::Scalar c = a.transpose()*(b-a);
        typename Kernel::Scalar bas = (b-a).squaredNorm(), as = a.squaredNorm();

        if(c<0) {
            beta = -c+std::sqrt(std::pow(c,2)+bas*(std::pow(delta,2)-as));
            beta /= bas;
        }
        else {
            beta = std::pow(delta,2)-as;
            beta /= c+std::sqrt(std::pow(c,2) + bas*(std::pow(delta,2)-as));
        };

        // and update h_dl and dL with beta
        h_dl = alpha*h_sd + beta*(b-a);

#ifdef USE_LOGGING
        BOOST_LOG_SEV(log, iteration) << "Dogleg step"; 
        if(!boost::math::isfinite(c)) {
            BOOST_LOG_SEV(log, error) << "Unnormal dogleg c detected: "<<c;
        }

        if(!boost::math::isfinite(bas)) {
            BOOST_LOG_SEV(log, error) << "Unnormal dogleg bas detected: "<<bas;
        }

        if(!boost::math::isfinite(beta)) {
            BOOST_LOG_SEV(log, error) << "Unnormal dogleg beta detected: "<<beta;
        }

#endif
    }
};

template<typename Kernel>
void Dogleg<Kernel>::execute() {

    clock_t start = clock();

    /*
    if(!m_system.isValid())
        throw solving_error() <<  boost::errinfo_errno(5) << error_message("invalid equation system");
    */
    
    F_old.resize(m_system.equationCount());
    g.resize(m_system.equationCount());
    J_old.resize(m_system.equationCount(), m_system.parameterCount());

    m_equations->execute();
    
#ifdef USE_LOGGING
    BOOST_LOG_SEV(log, solving) << "initial jacobi: "<<std::endl<<m_system.jacobi()<<std::endl
                                << "residual: "<<m_system.residuals().transpose()<<std::endl
                                << "maximal differential: "<<m_system.jacobi().template lpNorm<Eigen::Infinity>();
#endif
    
    //sys.removeLocalGradientZeros(true);
    m_equations->execute();
    //sys.removeLocalGradientZeros(false);
    
#ifdef USE_LOGGING
    BOOST_LOG_SEV(log, solving) << "LGZ jacobi: "<<std::endl<<m_system.jacobi()<<std::endl
                                << "maximal differential: "<<m_system.jacobi().template lpNorm<Eigen::Infinity>();
#endif

    err = m_system.residuals().norm();

    F_old = m_system.residuals();
    J_old = m_system.jacobi();

    g = m_system.jacobi().transpose()*(m_system.residuals());

    // get the infinity norm fx_inf and g_inf
    g_inf = g.template lpNorm<Eigen::Infinity>();
    fx_inf = m_system.residuals().template lpNorm<Eigen::Infinity>();

    delta=5;
    nu=2.;
    iter=0;
    stop=0;
    reduce=0;
    unused=0;
    counter=0;

    int maxIterNumber = 1000;
    typename Kernel::Scalar pr = 1e-6;//*m_system.Scaling;
    typename Kernel::Scalar diverging_lim = 1e6*err + 1e12;

    do {

        // check if finished
        if(fx_inf <= pr)  // Success
            stop = 1;
        else if(g_inf <= tolg*pr)
            throw solving_error() <<  boost::errinfo_errno(2) << error_message("g infinity norm smaller below limit");
        else if(delta <= tolx*pr)
            throw solving_error() <<  boost::errinfo_errno(3) << error_message("step size below limit");
        else if(iter >= maxIterNumber)
            throw solving_error() <<  boost::errinfo_errno(4) << error_message("maximal iterations reached");
        else if(!boost::math::isfinite(err))
            throw solving_error() <<  boost::errinfo_errno(5) << error_message("error is inf or nan");
        else if(err > diverging_lim)
            throw solving_error() <<  boost::errinfo_errno(6) << error_message("error diverged");


        // see if we are already finished
        if(stop)
            break;

        typename Kernel::Scalar err_new;
        typename Kernel::Scalar dF=0, dL=0;
        typename Kernel::Scalar rho;

        //get the update step
        calculateStep(g, m_system.jacobi(), m_system.residuals(), h_dl, delta);

#ifdef USE_LOGGING
        BOOST_LOG_SEV(log, iteration) << "Step in iter "<<iter<<std::endl
                                      << "Step: "<<h_dl.transpose()<<std::endl
                                      << "Jacobi: "<<m_system.jacobi()<<std::endl
                                      << "Residual: "<<m_system.residuals().transpose();
#endif

        // calculate the linear model
        dL = m_system.residuals().norm() - (m_system.residuals() + m_system.jacobi()*h_dl).norm();

        // get the new values
        m_system.parameter() += h_dl;
        m_equations->execute();

#ifdef USE_LOGGING

        if(!boost::math::isfinite(m_system.residuals().norm())) {
            BOOST_LOG_SEV(log, error) << "Unnormal residual detected: "<<m_system.residuals().norm();
        }

        if(!boost::math::isfinite(m_system.jacobi().sum())) {
            BOOST_LOG_SEV(log, error) << "Unnormal jacobi detected: "<<m_system.jacobi().sum();
        }

#endif

        //calculate the translation update ratio
        err_new = m_system.residuals().norm();
        dF = err - err_new;
        rho = dF/dL;

        if(dF<=0 || dL<=0)
            rho = -1;

        // update delta
        if(rho>0.85) {
            delta = std::max(delta,2*h_dl.norm());
            nu = 2;
        }
        else if(rho < 0.25) {
            delta = delta/nu;
            nu = 2*nu;
        }

#ifdef USE_LOGGING
        BOOST_LOG_SEV(log, iteration)<<"Result of step dF: "<<dF<<", dL: "<<dL<<std::endl
                                     << "New Residual: "<< m_system.residuals().transpose()<<std::endl;
#endif

        if(dF > 0 && dL > 0) {
/*
            //see if we got too high differentials
            if(m_system.jacobi().template lpNorm<Eigen::Infinity>() > 2) {
#ifdef USE_LOGGING
                BOOST_LOG_SEV(log, iteration)<< "High differential detected: "<<m_system.jacobi().template lpNorm<Eigen::Infinity>()<<" in iteration: "<<iter;
#endif
                rescale();
                m_equations->execute();
            }
            //it can also happen that the differentials get too small, however, we cant check for that
            else if(iter>1 && (counter>50)) {
                rescale();
                m_equations->execute();
                counter = 0;
            }
*/
            F_old = m_system.residuals();
            J_old = m_system.jacobi();

            err = m_system.residuals().norm();
            g = m_system.jacobi().transpose()*(m_system.residuals());

            // get infinity norms
            g_inf = g.template lpNorm<Eigen::Infinity>();
            fx_inf = m_system.residuals().template lpNorm<Eigen::Infinity>();
        }
        else {
#ifdef USE_LOGGING
            BOOST_LOG_SEV(log, iteration)<< "Reject step in iter "<<iter<<", dF: "<<dF<<", dL: "<<dL;
#endif
            m_system.residuals() = F_old;
            m_system.jacobi() = J_old;
            m_system.parameter() -= h_dl;
            unused++;
        }

        iter++;
        counter++;
    }
    while(!stop);


    clock_t end = clock();
    time = (double(end-start) * 1000.) / double(CLOCKS_PER_SEC);

#ifdef USE_LOGGING
    BOOST_LOG_SEV(log, solving)<<"Done solving: "<<err<<", iter: "<<iter<<", unused: "<<unused<<", reason:"<< stop;
    BOOST_LOG_SEV(log, solving)<< "final jacobi: "<<std::endl<<m_system.jacobi()
                               << "residual: "<<m_system.residuals().transpose()<<std::endl
                               << "maximal differential: "<<m_system.jacobi().template lpNorm<Eigen::Infinity>();
#endif
}

};

struct DummyKernel : public numeric::KernelBase {

    typedef int  Scalar;
};

template<typename NumericType, template<class> class Nonlinear = numeric::Dogleg>
struct Eigen3Kernel : public numeric::KernelBase {

    //the number type we use throughout the system
    typedef NumericType   Scalar;


private:
    //Nonlinear< Eigen3Kernel<Scalar, Nonlinear> > m_solver;

};




}//dcm


#endif //GCM_KERNEL_H






