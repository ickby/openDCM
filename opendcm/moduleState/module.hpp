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

#define DCM_USE_PARSER

#include <iosfwd>
#include <boost/shared_ptr.hpp>

#include <opendcm/core/object.hpp>
#include <opendcm/core/property.hpp>
#include <opendcm/core/clustergraph.hpp>
#include <opendcm/modulePart/module.hpp>
#include "traits.hpp"
#include "indent.hpp"
#include "generator.hpp"

#include <boost/mpl/int.hpp>
#include <boost/mpl/for_each.hpp>
#include <boost/fusion/container/vector/convert.hpp>
#include <boost/fusion/include/as_vector.hpp>
#include <boost/fusion/include/at.hpp>
#include <boost/fusion/include/at_c.hpp>

#include <boost/spirit/include/karma.hpp>
#include <boost/spirit/include/karma_rule.hpp>
#include <boost/spirit/include/karma_grammar.hpp>
#include <boost/spirit/include/karma_generate.hpp>
#include <boost/spirit/include/phoenix.hpp>

namespace karma = boost::spirit::karma;
namespace phx = boost::phoenix;
namespace fusion = boost::fusion;

namespace dcm {

struct ModuleState {

    template<typename Sys>
    struct type {

	typedef Unspecified_Identifier Identifier;

        template<typename Derived, typename Sig>
        void generate_properties(std::ostream stream, boost::shared_ptr< dcm::Object<Sys, Derived, Sig> > obj)  {

            typedef typename dcm::Object<Sys, Derived, Sig>::Sequence Properties;
        };

        struct inheriter {


            inheriter() :  m_this((Sys*) this) {}
            Sys* m_this;

            void saveState(std::ostream& stream) {

                typedef std::ostream_iterator<char> iterator_type;
                typedef typename boost::graph_traits<typename Sys::Cluster>::vertex_iterator viter;

                boost::iostreams::filtering_ostream indent_stream;
                indent_stream.push(indent_filter());
                indent_stream.push(stream);

                std::ostream_iterator<char> out(indent_stream);
                generator<iterator_type, Sys> gen(*m_this);

                karma::generate(out, gen, m_this->m_cluster);
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





