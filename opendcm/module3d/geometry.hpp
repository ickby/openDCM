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

#ifndef DCM_GEOMETRY_3D_H
#define DCM_GEOMETRY_3D_H

#include <opendcm/core/geometry.hpp>

namespace dcm {
namespace tag {

struct point3D  {
    typedef mpl::int_<3>  parameters;
    typedef mpl::int_<1>  transformations;
    typedef weight::point weight; 
};

struct line3D  {
    typedef mpl::int_<6> parameters;
    typedef mpl::int_<2> transformations;
    typedef weight::line weight; 
};

struct plane  {
    typedef mpl::int_<6>  parameters;
    typedef mpl::int_<2>  transformations;
    typedef weight::plane weight; 
};
}

namespace modell {
  
  struct XYZ {
    /*Modell XYZ: 
     * 0 = X;
     * 1 = Y;
     * 2 = Z;
     */    
    template<typename Scalar, typename Accessor, typename Vector, typename Type>
    void extract(Type& t, Vector& v) {
      Accessor a;
      v(0) = a.get<Scalar, 0>(t);
      v(1) = a.get<Scalar, 1>(t);
      v(2) = a.get<Scalar, 2>(t);
    }
  };
  
}

}

#endif //DCM_GEOMETRY_3D_H
