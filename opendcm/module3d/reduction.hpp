/*
    openDCM, dimensional constraint manager
    Copyright (C) 2016  Stefan Troeger <stefantroeger@gmx.net>

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

#ifndef DCM_MODULE3D_REDUCTION_H
#define DCM_MODULE3D_REDUCTION_H

#include "opendcm/core/reduction.hpp"
#include "geometry.hpp"

namespace dcm {
namespace module3d {
    
template<typename  Kernel>
struct FixedPoint : public numeric::DependendGeometry<Kernel, geometry::Point3, geometry::Point3> {
     
    using Inherited = numeric::DependendGeometry<Kernel, geometry::Point3, geometry::Point3>;

    CALCULATE() {
        Inherited::calculate();
        Inherited::output() = Inherited::input();
    };
};

template<typename Final>
void setupPointPointReduction(reduction::EdgeReductionGraph* graph) {

    typedef typename Final::Kernel Kernel;   
    typedef reduction::ConstraintEqualValue<Final, dcm::Distance, 0> ZeroDistance;
    
    auto fixPoint = graph->template getTreeNode<FixedPoint<Kernel>>();
    graph->sourceNode()->template connectConditional<ZeroDistance>(fixPoint, [](const dcm::Distance& dist){});      
};

}//module3d
}//dcm

#endif //DCM_MODULE3D_REDUCTION_H