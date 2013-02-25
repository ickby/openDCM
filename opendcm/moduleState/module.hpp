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

#ifndef DCM_MODULE_STATE_H
#define DCM_MODULE_STATE_H

#include <iosfwd>
#include "indent.hpp"

//forward declare the generate function so that we don't need to parse the spirit header when we externalize
//and when we are in a compilation unit which does not compile the generator. Spirit header consume LOTS of
//memory
#ifdef USE_EXTERNAL
namespace dcm {
template<typename Sys>
void generate(Sys* m_this, std::ostream& stream);
}
#else
#include "generator.hpp"
#endif

namespace dcm {

struct ModuleState {

    template<typename Sys>
    struct type {

        typedef Unspecified_Identifier Identifier;

        struct inheriter {

            inheriter() :  m_this((Sys*) this) {}
            Sys* m_this;

            void saveState(std::ostream& stream) {
                generate(m_this, stream);
            };

            void loadState(std::istream& stream) {

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

#endif //DCM_MODULE_STATE_H





