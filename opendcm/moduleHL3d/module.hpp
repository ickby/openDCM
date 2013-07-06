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

#ifndef GCM_MODULE_HLG3D_H
#define GCM_MODULE_HLG3D_H

#include <opendcm/core.hpp>
#include <opendcm/module3d.hpp>

#include "defines.hpp"
#include "geometry.hpp"
#include "generator.hpp"

#include <boost/mpl/if.hpp>
#include <boost/mpl/map.hpp>
#include <boost/type_traits.hpp>

#include <boost/preprocessor.hpp>
#include <boost/preprocessor/repetition/repeat.hpp>
#include <boost/preprocessor/cat.hpp>
#include <boost/preprocessor/repetition/enum.hpp>
#include <boost/preprocessor/repetition/enum_params.hpp>
#include <boost/preprocessor/repetition/enum_trailing_params.hpp>
#include <boost/preprocessor/repetition/enum_binary_params.hpp>


#define APPEND_SINGLE(z, n, data) \
    g_ptr = details::converter<BOOST_PP_CAT(Arg,n), Geometry3D, mpl::true_>::apply(BOOST_PP_CAT(arg,n), m_this); \
    if(!g_ptr) { \
      hlg_ptr = details::converter<BOOST_PP_CAT(Arg,n), HLGeometry3D, mpl::false_>::apply(BOOST_PP_CAT(arg,n), m_this); \
      if(!hlg_ptr) \
	throw creation_error() <<  boost::errinfo_errno(216) << error_message("could not handle input"); \
      else \
	data->append(hlg_ptr);\
    } \
    else { \
      data->append(g_ptr); \
    } \
 

#define CREATE_DEF(z, n, data) \
    template < \
    typename Generator  \
    BOOST_PP_ENUM_TRAILING_PARAMS(n, typename Arg) \
    > \
    boost::shared_ptr<HLGeometry3D> createHLGeometry3D( \
                     BOOST_PP_ENUM_BINARY_PARAMS(n, Arg, const& arg) \
                   );

#define CREATE_DEC(z, n, data) \
    template<typename ID> \
    template<typename Sys> \
    template <\
    typename Generator \
    BOOST_PP_ENUM_TRAILING_PARAMS(n, typename Arg)\
    > \
    boost::shared_ptr<typename ModuleHL3D<ID>::template type<Sys>::HLGeometry3D> \
    ModuleHL3D<ID>::type<Sys>::inheriter_base::createHLGeometry3D( \
            BOOST_PP_ENUM_BINARY_PARAMS(n, Arg, const& arg) \
                                              ) \
    { \
      typedef typename system_traits<Sys>::template getModule<details::m3d>::type module3d; \
      typedef typename module3d::template type<Sys>::Geometry3D Geometry3D; \
      boost::shared_ptr<Geometry3D> g_ptr; \
      boost::shared_ptr<HLGeometry3D> hlg_ptr; \
      boost::shared_ptr<HLGeometry3D> ptr = boost::shared_ptr<HLGeometry3D>(new HLGeometry3D(*m_this)); \
      BOOST_PP_REPEAT(n, APPEND_SINGLE, ptr) \
      ptr->template init<Generator>();\
      return ptr;\
    };


namespace mpl = boost::mpl;

namespace dcm {

namespace details {

//return always a geometry3d pointer struct, no matter whats the supplied type
template<typename T, typename R, typename Create>
struct converter {
    template<typename Sys>
    static boost::shared_ptr<R> apply(T const& t, Sys* sys) {
      return boost::shared_ptr<R>();
    };
};
template<typename T, typename R>
struct converter< T, R, mpl::true_> {
    template<typename Sys>
    static boost::shared_ptr<R> apply(T const& t, Sys* sys) {
      return sys->createGeometry3D(t);
    };
};
template<typename R, typename Create>
struct converter< boost::shared_ptr<R>, R, Create > {
    template<typename Sys>
    static boost::shared_ptr<R> apply(boost::shared_ptr<R> t, Sys* sys) {
        return t;
    };
};
template<typename T, typename R, typename Create>
struct converter< boost::shared_ptr<T>, R, Create >  {
    template<typename Sys>
    static boost::shared_ptr<R> apply(boost::shared_ptr<T> t, Sys* sys) {
        return  boost::shared_ptr<R>();
    };
};

}//details


template<typename ID = No_Identifier>
struct ModuleHL3D {

    template<typename Sys>
    struct type : details::mhl3d {

        //forward declare
        struct inheriter_base;

        struct HLGeometry3D : public Object<Sys, HLGeometry3D, mpl::map0<> > {

            HLGeometry3D(Sys& sys) : Object<Sys, HLGeometry3D, mpl::map0<> >(sys) {};

            template<typename tag>
            typename details::tag_traits<tag>::return_type get() {

            };

        protected:

            boost::shared_ptr< details::HLGeneratorBase<Sys> > m_generator;

            //traits are only accessible in subclass scope
            BOOST_MPL_ASSERT((typename system_traits<Sys>::template getModule<details::m3d>::has_module));
            typedef typename system_traits<Sys>::template getModule<details::m3d>::type module3d;
            typedef typename module3d::template type<Sys>::Geometry3D Geometry3D;
            typedef typename module3d::template type<Sys>::Constraint3D Constraint3D;


            using Object<Sys, HLGeometry3D, mpl::map0<> >::m_system;
            typedef Object<Sys, HLGeometry3D, mpl::map0<> > base;

            std::vector<boost::shared_ptr<Geometry3D> > m_geometrys;
            std::vector<boost::shared_ptr<HLGeometry3D> > m_hlgs;
            std::vector<boost::shared_ptr<Constraint3D> > m_constraints;

            template<typename generator>
            void init() {
                m_generator = boost::shared_ptr<details::HLGeneratorBase<Sys> >(new typename generator::template type<Sys>);
                m_generator->set(&m_geometrys, &m_hlgs, &m_constraints);

                if(!m_generator->check())
                    throw creation_error() <<  boost::errinfo_errno(210) << error_message("not all needd types for high level geometry present");

                m_generator->init();
            };


            boost::shared_ptr<HLGeometry3D> append(boost::shared_ptr<Geometry3D> g) {
                m_geometrys.push_back(g);
                return base::shared_from_this();
            };
            boost::shared_ptr<HLGeometry3D> append(boost::shared_ptr<HLGeometry3D> g) {
                m_hlgs.push_back(g);
                return base::shared_from_this();
            };

            friend struct inheriter_base;
	    friend struct Object<Sys, HLGeometry3D, mpl::map0<> >;
        };

        //inheriter for own functions
        struct inheriter_base {

            inheriter_base() : m_this((Sys*)this) {};

        private:
            Sys* m_this;

        public:
            //with no vararg templates before c++11 we need preprocessor to create the overloads of create we need
            BOOST_PP_REPEAT(5, CREATE_DEF, ~)

        };
        struct inheriter_id : public inheriter_base {};
        struct inheriter : public mpl::if_<boost::is_same<Identifier, No_Identifier>, inheriter_base, inheriter_id>::type {};

        //needed typedefs
        typedef ID Identifier;
        typedef mpl::vector0<> properties;
        typedef mpl::vector0<> objects;

        //needed static functions
        static void system_init(Sys& sys) {};
        static void system_copy(const Sys& from, Sys& into) {};

    };
};

/*****************************************************************************************************************/
/*****************************************************************************************************************/
/*****************************************************************************************************************/
/*****************************************************************************************************************/

BOOST_PP_REPEAT(5, CREATE_DEC, ~)


}//dcm

#endif //GCM_MODULE_HLG3D_H
