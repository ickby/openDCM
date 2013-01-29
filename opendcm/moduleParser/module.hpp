/*
    openDCM, dimensional constraint manager
    Copyright (C) 2013  Stefan Troeger <stefantroeger@gmx.net>

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

#ifndef DCM_MODULE_PARSER_H
#define DCM_MODULE_PARSER_H

#define DCM_USE_PARSER

#include <iosfwd>
#include <boost/shared_ptr.hpp>
#include <opendcm/core/object.hpp>


namespace dcm {

struct ModuleParser {

    template<typename Sys>
    struct type {

        template<typename Derived, typename Sig>
        struct propperty_generator {

            std::ostream& stream;
            boost::shared_ptr< dcm::Object<Sys, Derived, Sig> > object;

            propperty_generator(std::ostream& s,
                                boost::shared_ptr< dcm::Object<Sys, Derived, Sig> > o) : stream(s), object(o) {};

            template< typename T >
            void operator()(T x) {
	      
            };
        };

        template<typename Derived, typename Sig>
        void generate_properties(std::ostream stream, boost::shared_ptr< dcm::Object<Sys, Derived, Sig> > obj)  {

            typedef typename dcm::Object<Sys, Derived, Sig>::Sequence Properties;
        };

        struct inheriter {

            inheriter() :  m_this((Sys*) this) {}
            Sys* m_this;

            void saveTo(std::ostream stream) {

            };
        };


        //nothing to add to the objects and properties of other modules
        typedef mpl::vector0<>  properties;
        typedef mpl::vector0<>  objects;

        //nothing to do on startup
        static void system_init(Sys& sys) {};
    };
};

}

#endif //DCM_MODULEPARSER_H




