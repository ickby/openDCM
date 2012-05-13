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
#include <eigen3/Eigen/Geometry>

namespace dcm {
  
  using namespace Eigen;

template<typename Nt>
class Kernel {

  double precision;
  
public:
    typedef Nt number_type;
    typedef Matrix<Nt, 3, 1> Vector3;
    typedef Matrix<Nt, 3, 3> Matrix3;
    
    Kernel() : precision(0.001) {};

    template <typename DerivedA,typename DerivedB>
    static typename DerivedA::Scalar isSame(const MatrixBase<DerivedA>& p1,const MatrixBase<DerivedB>& p2) {
        return ((p1-p2).squaredNorm() < 0.001);
    }
    template <typename DerivedA,typename DerivedB>
    static typename DerivedA::Scalar isOpposite(const MatrixBase<DerivedA>& p1,const MatrixBase<DerivedB>& p2) {
        return ((p1+p2).squaredNorm() < 0.001);
    }
};

}

#endif //DCM_KERNEL_H
