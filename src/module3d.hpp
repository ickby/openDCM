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

#ifndef DCM_GEOMETRY3D_H
#define DCM_GEOMETRY3D_H

#include <boost/mpl/vector.hpp>
#include <boost/mpl/less.hpp>
#include <boost/mpl/if.hpp>
#include <boost/mpl/size.hpp>

#include <boost/static_assert.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/preprocessor/iteration/local.hpp>
#include <boost/variant.hpp>

#include "object.hpp"
#include "geometry.hpp"
#include "constraint3d.hpp"

namespace mpl = boost::mpl;

namespace dcm {

namespace details {

template<typename seq, typename t>
struct distance {
    typedef typename mpl::find<seq, t>::type iterator;
    typedef typename mpl::distance<typename mpl::begin<seq>::type, iterator>::type type;
    BOOST_MPL_ASSERT(( mpl::not_< boost::is_same<iterator, typename mpl::end<seq>::type > > ));
};

}

template<typename Typelist>
struct Module3D {

    template<typename Sys>
    struct type {
        class Constraint3D;

        class Geometry3D : public Object<Sys, Geometry3D, mpl::map<> > {
            typedef typename boost::make_variant_over< Typelist >::type Variant;

        public:
            template<typename T>
            Geometry3D(T geometry, Sys& system) : Object<Sys, Geometry3D, mpl::map<> >(system),
                    m_geometry(geometry) {

                m_type = details::distance<Typelist, T>::type::value;
            };

        protected:
            Variant m_geometry;
            Storage m_storage;
            int m_type;

            friend class Constraint3D;
        };

        typedef boost::shared_ptr<Geometry3D> Geom;

        //type erasure container for constraints
        class Constraint3D : public Object<Sys, Constraint3D, mpl::map<> > {

        public:

            Constraint3D(Sys& system, Geom f, Geom s) : Object<Sys, Constraint3D, mpl::map<> >(system),
                    first(f), second(s), content(0) {  };

            ~Constraint3D()  {
                delete content;
            }

            template<template<typename, typename> class T>
            void setType() {
                creator<T> creator;
                boost::apply_visitor(creator, first->m_geometry, second->m_geometry);
                content = creator.p;
            };

            double calculate() {
                return content->calculate(first->m_storage, second->m_storage);
            };

        protected:

            struct placeholder  {

                virtual ~placeholder() {}
                virtual double calculate(Storage&, Storage&) const = 0;
                virtual placeholder* resetConstraint(Geom first, Geom second) const = 0;
            };

            template< template<typename,typename> class T1, typename T2, typename T3>
            struct holder : public placeholder  {

                holder(const T1<T2,T3> & value)
                        : held(value)   {}

                virtual double calculate(Storage& f, Storage& s) const {
                    return held.calculate(f,s);
                };
                T1<T2,T3>  held;

                virtual placeholder* resetConstraint(Geom first, Geom second) const {
                    creator<T1> creator;
                    boost::apply_visitor(creator, first->m_geometry, second->m_geometry);
                    return creator.p;
                };
            };


            template<template<typename, typename> class T>
            struct creator : public boost::static_visitor<void> {

                template<typename T1, typename T2>
                void operator() (const T1 &, const T2 &) {
                    typedef T<typename geometry_traits<T1>::tag, typename geometry_traits<T2>::tag> type;
                    p = new holder< T, typename geometry_traits<T1>::tag, typename geometry_traits<T2>::tag >( type() );
                };
                placeholder* p;
            };

            placeholder * content;
            Geom first, second;

        };

        typedef boost::shared_ptr<Constraint3D> Cons;

        typedef mpl::vector<Geometry3D, Constraint3D> objects;

        struct inheriter {
            template<typename T>
            Geom createGeometry3D(T geom ) {
                return Geom(new Geometry3D(geom, *((Sys*)this)));
            };

            template<template<typename, typename> class T>
            Cons createConstraint3D(Geom first, Geom second) {
                Cons c(new Constraint3D(*((Sys*)this), first, second));
                c->template setType<T>();
                return c;
            };
        };
        typedef mpl::vector<>  	properties;

    };
};

}

#endif //DCM_GEOMETRY3D_H
