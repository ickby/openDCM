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

#ifndef DCM_MODULE_3D_H
#define DCM_MODULE_3D_H

#include "opendcm/core/object.hpp"
#include "opendcm/core/module.hpp"
#include "opendcm/core/utilities.hpp"

namespace dcm {
    
template<typename ... types>
struct Module3D {

    typedef boost::mpl::int_<5> ID;

    template<typename Final, typename Stacked>
    struct type : public Stacked {

        struct Geometry3D : public Stacked::ObjectBase, utilities::Variant<types...> {

            Geometry3D() : Stacked::ObjectBase( Final::template objectTypeID<typename Final::Geometry3D>::ID::value ) {};
            
            //DCM_OBJECT_ADD_PROPERTIES( Final, (TestProperty1)(TestProperty2) )
        };
        
        DCM_MODULE_ADD_OBJECTS(Stacked, (Geometry3D))
    };
    
};

}//dcm

#include "imp/module_imp.hpp"

#endif //DCM_GEOMETRY3D_H







