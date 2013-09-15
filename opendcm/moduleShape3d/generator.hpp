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

#ifndef GCM_GENERATOR_SHAPE3D_H
#define GCM_GENERATOR_SHAPE3D_H

#include <opendcm/core/defines.hpp>
#include <opendcm/core/geometry.hpp>
#include <opendcm/module3d/defines.hpp>
#include "geometry.hpp"
#include "fixed.hpp"

#include <boost/exception/errinfo_errno.hpp>

namespace dcm {

namespace details {

template<typename Sys>
struct HLGeneratorBase {

    BOOST_MPL_ASSERT((typename system_traits<Sys>::template getModule<details::m3d>::has_module));
    typedef typename system_traits<Sys>::template getModule<details::m3d>::type module3d;
    typedef typename module3d::Geometry3D Geometry3D;
    typedef typename module3d::Constraint3D Constraint3D;
    typedef typename system_traits<Sys>::template getModule<details::mshape3d>::type modulehl3d;
    typedef typename modulehl3d::Shape3D Shape3D;

    Sys* m_system;
    boost::shared_ptr<Shape3D> m_shape;
    std::vector<boost::shared_ptr<Geometry3D> >*   m_geometrys;
    std::vector<boost::shared_ptr<Shape3D> >* m_hlgs;
    std::vector<boost::shared_ptr<Constraint3D> >* m_constraints;

    HLGeneratorBase(Sys* system) : m_system(system) {};
    virtual ~HLGeneratorBase() {};

    void set(boost::shared_ptr<Shape3D> shape,
             std::vector<boost::shared_ptr<Geometry3D> >*   geometrys,
             std::vector<boost::shared_ptr<Shape3D> >* hlgs,
             std::vector<boost::shared_ptr<Constraint3D> >* constraints) {

        m_shape = shape;
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
    virtual boost::shared_ptr<Shape3D> getOrCreateHLG3d(int type) = 0;
};

} //details


struct dummy_generator {

    template<typename Sys>
    struct type : public details::HLGeneratorBase<Sys> {

        type(Sys* system) : details::HLGeneratorBase<Sys>(system) {};

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
        virtual boost::shared_ptr<typename details::HLGeneratorBase<Sys>::Shape3D> getOrCreateHLG3d(int type) {
            throw creation_error() <<  boost::errinfo_errno(213) << error_message("dummy generator has no high level geometry to access");
        };
    };
};

//test generator
struct segment3D {

    template<typename Sys>
    struct type : public dcm::details::HLGeneratorBase<Sys> {

        typedef dcm::details::HLGeneratorBase<Sys> base;
        typedef typename Sys::Kernel Kernel;
        using typename base::Geometry3D;
	using typename base::Constraint3D;

        type(Sys* system) : details::HLGeneratorBase<Sys>(system) {};

        //check if all needed parts are supplied, a segment needs 2 points
        virtual bool check() {

            //even we have a real geometry segment
            if(base::m_shape->getGeneralType() == tag::weight::segment::value)
                return true;

            //or two point geometries
            if(base::m_geometrys->size() == 2)
                return true;

            return false;
        };
        //initialise all relations between the geometrys
        virtual void init() {

            if(base::m_shape->getGeneralType() == dcm::tag::weight::segment::value) {
                //we have a segment, lets link the two points to it


            }

            if(base::m_geometrys->size() == 2) {
                //we have two points, lets get them
                boost::shared_ptr<Geometry3D> g1 = base::m_geometrys->operator[](0);
                boost::shared_ptr<Geometry3D> g2 = base::m_geometrys->operator[](1);

                //possibility 1: two points. we add a segment line an link the point in
                if(g1->getGeneralType() == tag::weight::point::value || g2->getGeneralType() == tag::weight::point::value) {

                    //construct our segment value
                    typename Kernel::Vector val(6);
                    val.head(3) = g1->getValue();
                    val.tail(3) = g2->getValue();

                    //and create a segment geometry we use as line
                    boost::shared_ptr<Geometry3D> g3 = base::m_system->createGeometry3D();
                    g3->template setValue<tag::segment3D>(val);

                    //link the points to our new segment
                    g1->template linkTo<tag::point3D>(g3, 0);
		    g2->template linkTo<tag::point3D>(g3, 3);

                    //add the fix constraints to show our relation
		    boost::shared_ptr<Constraint3D> c1 = base::m_system->createConstraint3D(g1,g3, details::fixed);
		    boost::shared_ptr<Constraint3D> c2 = base::m_system->createConstraint3D(g1,g3, details::fixed);
		    c1->disable(); //required by fixed constraint
		    c2->disable(); //requiered by fixed constraint

                }
                else
                    throw creation_error() <<  boost::errinfo_errno(501) << error_message("Wrong geometries for segment construction");
            };
        };
        //get geometry3d for optional types (e.g. midpoints)
        virtual boost::shared_ptr<typename dcm::details::HLGeneratorBase<Sys>::Geometry3D> getOrCreateG3d(int type) {
            return boost::shared_ptr<typename dcm::details::HLGeneratorBase<Sys>::Geometry3D>();
        };
        //get hlgeometry3d for optional types
        virtual boost::shared_ptr<typename dcm::details::HLGeneratorBase<Sys>::Shape3D> getOrCreateHLG3d(int type) {
            return boost::shared_ptr<typename dcm::details::HLGeneratorBase<Sys>::Shape3D>();
        };
    };
};

} //dcm


#endif //GCM_GENERATOR_SHAPE3D_H
