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

#ifndef DCM_CONSTRAINT3D_H
#define DCM_CONSTRAINT3D_H

namespace dcm {

template<typename T1, typename T2>
struct Coincident3D {

    //TODO: throw as someone called a unsupported geometry combination with this constraint
};

template<>
struct Coincident3D<tag_point3D, tag_point3D> {

    template<typename Kernel>
    typename Kernel::number_type calculate(typename Kernel::VectorMap& param1,
                                           typename Kernel::VectorMap& param2) const {};

    template<typename Kernel>
    typename Kernel::number_type calculateFirstFullGradient(typename Kernel::VectorMap& param1,
            typename Kernel::VectorMap& param2,
            typename Kernel::VectorMap& diffparam) const {};

    template<typename Kernel>
    typename Kernel::number_type calculateSecondFullGradient(typename Kernel::VectorMap& param1,
            typename Kernel::VectorMap& param2,
            typename Kernel::VectorMap& diffparam) const {};

    template<typename Kernel>
    void calculateFirstGradient(typename Kernel::VectorMap& param1,
                                typename Kernel::VectorMap& param2,
                                typename Kernel::VectorMap& grad) const {};

    template<typename Kernel>
    void calculateSecondGradient(typename Kernel::VectorMap& param1,
                                 typename Kernel::VectorMap& param2,
                                 typename Kernel::VectorMap& grad) const {};

};

}

#endif //DCM_CONSTRAINT3D_H
