/*
    openDCM, dimensional constraint manager
    Copyright (C) 2013  Stefan Troeger <stefantroeger@gmx.net>

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

#ifndef GCM_GENERATOR_HLG3D_H
#define GCM_GENERATOR_HLG3D_H

#include <opendcm/core/defines.hpp>
#include <boost/exception/errinfo_errno.hpp>

namespace dcm {

namespace details {

template<typename Sys>
struct HLGeneratorBase {

    BOOST_MPL_ASSERT((typename system_traits<Sys>::template getModule<details::m3d>::has_module));
    typedef typename system_traits<Sys>::template getModule<details::m3d>::type module3d;
    typedef typename module3d::template type<Sys>::Geometry3D Geometry3D;
    typedef typename module3d::template type<Sys>::Constraint3D Constraint3D;
    typedef typename system_traits<Sys>::template getModule<details::mhl3d>::type modulehl3d;
    typedef typename modulehl3d::template type<Sys>::HLGeometry3D HLGeometry3D;

    std::vector<boost::shared_ptr<Geometry3D> >*   m_geometrys;
    std::vector<boost::shared_ptr<HLGeometry3D> >* m_hlgs;
    std::vector<boost::shared_ptr<Constraint3D> >* m_constraints;

    HLGeneratorBase() {};
    virtual ~HLGeneratorBase() {};

    void set(std::vector<boost::shared_ptr<Geometry3D> >*   geometrys,
             std::vector<boost::shared_ptr<HLGeometry3D> >* hlgs,
             std::vector<boost::shared_ptr<Constraint3D> >* constraints) {

        m_geometrys = geometrys;
        m_hlgs = hlgs;
        m_constraints = constraints;
    };

    //check if all needed parts are supplied
    virtual bool check() = 0;
    //initialise all relations between the geometrys
    virtual void init() = 0;
    //get geometry3d for optional types (e.g. midpoints)
    virtual boost::shared_ptr<Geometry3D> getOrCreateG3d(int type) = 0;
    //get hlgeometry3d for optional types
    virtual boost::shared_ptr<HLGeometry3D> getOrCreateHLG3d(int type) = 0;
};

} //details


struct dummy_generator {

    template<typename Sys>
    struct type : public details::HLGeneratorBase<Sys> {

        //check if all needed parts are supplied
        virtual bool check() {
            throw creation_error() <<  boost::errinfo_errno(210) << error_message("not all needd types for high level geometry present");
        };
        //initialise all relations between the geometrys, throw on error
        virtual void init() {
            throw creation_error() <<  boost::errinfo_errno(211) << error_message("dummy generator can't create high level geometry");
        };
        //get geometry3d for optional types (e.g. midpoints)
        virtual boost::shared_ptr<typename details::HLGeneratorBase<Sys>::Geometry3D> getOrCreateG3d(int type) {
            throw creation_error() <<  boost::errinfo_errno(212) << error_message("dummy generator has no geometry to access");
        };
        //get hlgeometry3d for optional types
        virtual boost::shared_ptr<typename details::HLGeneratorBase<Sys>::HLGeometry3D> getOrCreateHLG3d(int type) {
            throw creation_error() <<  boost::errinfo_errno(213) << error_message("dummy generator has no high level geometry to access");
        };
    };
};

} //dcm


#endif //GCM_GENERATOR_HLG3D_H
