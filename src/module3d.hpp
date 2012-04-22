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

#ifndef NS2_GEOMETRY3D_H
#define NS2_GEOMETRY3D_H

#include <boost/mpl/vector.hpp>

#include "object.hpp"

namespace mpl = boost::mpl;

namespace dcm {

//template<typename T>
struct Module3D {

    template<typename Sys>
    struct type {

        struct Geometry3D : public Object<Sys, Geometry3D> {};
	typedef mpl::vector<Geometry3D> objects;

        struct inheriter {
	  boost::shared_ptr<Geometry3D> createGeometry3D() {
	    return boost::shared_ptr<Geometry3D>(new Geometry3D);
	  };
	};
        typedef mpl::map<>	obj_properties;
        typedef mpl::vector<>  	graph_properties;

    };
};

}

#endif //NS2_GEOMETRY3D_H
