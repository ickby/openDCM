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
    GNU General Public License for more detemplate tails.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/

#ifndef GCM_COINCIDENT_H
#define GCM_COINCIDENT_H

#include "geometry.hpp"

namespace gcm {


template<typename Kernel, typename Tag1, typename Tag2>
struct Coincident3D {

    typedef typename Kernel::number_type Scalar;
    typedef typename Kernel::VectorMap   Vector;
    
    Scalar getEquationScaling(typename Kernel::Vector& local1, typename Kernel::Vector& local2) {
      assert(false);
    }

    Scalar calculate(Vector& param1,  Vector& param2)  {
        assert(false);
    };

    Scalar calculateGradientFirst(Vector& param1,  Vector& param2, Vector& dparam1) {
        assert(false);
    };

    Scalar calculateGradientSecond(Vector& param1,  Vector& param2, Vector& dparam2)  {
        assert(false);
    };

    void calculateGradientFirstComplete(Vector& param1,  Vector& param2, Vector& gradient) {
        assert(false);
    };

    void calculateGradientSecondComplete(Vector& param1,  Vector& param2, Vector& gradient) {
        assert(false);
    };

};

}

#endif //GCM_COINCIDENT_H