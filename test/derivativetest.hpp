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
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License along
    with this library; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/

#ifndef DCM_DERIVATIVE_TEST
#define DCM_DERIVATIVE_TEST

#include <boost/test/unit_test.hpp>
#include <opendcm/core/geometry.hpp>

#define DELTA  1e-6
#define MAXDIV 1e-3
#define STEP   0.1
#define RANGE  5

struct DerivativeTest {

    
    template<typename Kernel, template<class, bool> class Base>
    struct DerivativeCalculator {
        
        typedef Base<Kernel, false> Geometry;
        
        Geometry &left, &right, &diff;
        
        DerivativeCalculator(Geometry& l, Geometry& r, Geometry& d) : left(l), right(r), diff(d) {};
        
        template<typename T>
        void operator()(T& i) const {
            
            auto& ls = fusion::at<T>(left.m_storage);
            auto& rs = fusion::at<T>(right.m_storage);
            auto& d = fusion::at<T>(diff.m_storage);
            
            calc(ls, rs, d);
        };
        
        template<typename T>
        void calc(Eigen::MatrixBase<T>& v1, Eigen::MatrixBase<T>& v2,
            Eigen::MatrixBase<T>& d) const {
            
            Eigen::MatrixXd res = (v2-v1);
            d = res / (2*DELTA);
        };
        
        void calc(double* v1, double* v2, double* d) const {
            *d = (*v2 - *v1)  / (2*DELTA);
        };
    };
    
    template<typename Kernel, template<class, bool> class Base>
    struct Compare {
      
        Base<Kernel, false>& analytical;
        Base<Kernel, false>& numeric;
        int iteration, derivative;
        
        Compare(Base<Kernel, false>& a, Base<Kernel, false>& n, int i, int d) 
            : analytical(a), numeric(n), iteration(i), derivative(d) {};
        
        template<typename T>
        void operator()(T& i) {
            auto& ls = fusion::at<T>(analytical.m_storage);
            auto& rs = fusion::at<T>(numeric.m_storage);
            
            if(!calc(ls,rs)) {
                std::cout<<"storage "<< T::value << " in iteration " << iteration << " of derivative "<< derivative
                    <<std::endl<< "numeric: "<<std::endl
                    <<rs<<std::endl<<" vs analytic: "<<std::endl<<ls<<std::endl;
                throw(std::exception());
            }
        };
        
        template<typename T>
        bool calc(Eigen::MatrixBase<T>& v1, Eigen::MatrixBase<T>& v2) {
            return (v1-v2).norm() < MAXDIV;
        };
        
        bool calc(double* v1, double* v2) {
            return std::abs((*v1 - *v2)) < MAXDIV;
        };
    };
    
    template<typename Kernel, template<class, bool> class Base>
    static bool isCorrect(dcm::numeric::Geometry<Kernel, Base>& g, const std::function<void()>& recalc) {

        typedef typename dcm::numeric::Geometry<Kernel, Base>::Derivative Derivative;

        int d = 0;
        for(Derivative& der : g.derivatives()) {

            for(int i=0; i<int(RANGE/STEP); i++) {

                *(der.second.Value) += STEP - DELTA;
                recalc();
                //save the current value
                Base<Kernel, false> left;
                left.m_storage = g.m_storage;
                
                *(der.second.Value) += 2*DELTA;
                recalc();
                //save the current value
                Base<Kernel, false> right;
                right.m_storage = g.m_storage;

                //get the residual diff between the last two calculations
                Base<Kernel, false> diff;
                //mpl trickery to get a sequence counting from 0 to the size of stroage entries
                typedef mpl::range_c<int,0,
                        mpl::size<typename Base<Kernel, false>::StorageTypeSequence>::value> StorageRange;
                //now iterate that sequence so we can access all storage elements with knowing the position
                //we are at (that is important to access the correct derivative storage position too)
                mpl::for_each<StorageRange>( DerivativeCalculator<Kernel, Base>(left, right, diff) );
                
                //compare to the analytical diff
                Base<Kernel, false> analytic;
                analytic.m_storage = der.first.m_storage;
                mpl::for_each<StorageRange>(Compare<Kernel, Base>(analytic, diff, i, d));
                
            };
            ++d;
        };
        
        return true;
    };
};

#endif
