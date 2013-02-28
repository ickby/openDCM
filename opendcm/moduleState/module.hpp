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
#include "generator.hpp"
#include "parser.hpp"

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/support_istream_iterator.hpp>

namespace qi = boost::spirit::qi;

namespace dcm {

struct ModuleState {

    template<typename Sys>
    struct type {

        typedef Unspecified_Identifier Identifier;

        struct inheriter {

            inheriter() :  m_this((Sys*) this) {}
            Sys* m_this;

            void saveState(std::ostream& stream) {

                boost::iostreams::filtering_ostream indent_stream;
                indent_stream.push(indent_filter());
                indent_stream.push(stream);

                std::ostream_iterator<char> out(indent_stream);
                generator<Sys> gen(*m_this);

                karma::generate(out, gen, m_this->m_cluster);
            };

            void loadState(std::istream& stream) {
	      
                //disable skipping of whitespace
                stream.unsetf(std::ios::skipws);

                // wrap istream into iterator
                boost::spirit::istream_iterator begin(stream);
                boost::spirit::istream_iterator end;

                // use iterator to parse file data
		parser<Sys> par(*m_this);
		m_this->clear();
                qi::parse(begin, end, par, m_this->m_cluster);
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





