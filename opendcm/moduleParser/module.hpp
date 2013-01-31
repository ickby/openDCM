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
#include <opendcm/core/property.hpp>
#include "traits.hpp"
#include "container.hpp"
#include <opendcm/modulePart/module.hpp>

#include <boost/spirit/include/karma.hpp>
#include <boost/spirit/include/karma_rule.hpp>
#include <boost/spirit/include/karma_grammar.hpp>
#include <boost/spirit/include/karma_generate.hpp>

namespace karma = boost::spirit::karma;

namespace dcm {

struct ModuleParser {

    template<typename Sys>
    struct type {

        //generates the output for a specific property and pushes it into the supplied ostream.
        template< typename Property >
        struct propperty_generator {

            typedef std::back_insert_iterator<std::string> iterator_type;
            typedef typename dcm::parser_generator<Property, Sys, iterator_type>::generator generator;

            generator gen;
            std::ostream& stream;

            propperty_generator(std::ostream& s) : stream(s) {};

            void operator()(typename Property::type& t) {
                std::string s;
                std::back_insert_iterator<std::string> out(s);
                dcm::parser_generator<Property, Sys, iterator_type>::init(gen);
                boost::spirit::karma::generate(out, gen, t);
                stream << s;
            };
        };

        //generates the output for a specific property out of a property vector, the values of all
        //propertys must be supplied as fusion sequence of the property types. This allows the use
        //of mpl::for_each with this strcut as functor.
        template< typename PropertyVector >
        struct propperty_vector_generator {

            typedef std::back_insert_iterator<std::string> iterator_type;
            typedef typename details::pts<PropertyVector>::type fusion_sequence;

            fusion_sequence& sequence;
            std::ostream& stream;

            propperty_vector_generator(std::ostream& s, fusion_sequence& seq) : stream(s), sequence(seq) {};

            template<typename Property>
            void operator()(Property p) {

                typedef typename mpl::find<PropertyVector, Property>::type iterator;
                typedef typename mpl::distance<typename mpl::begin<PropertyVector>::type, iterator>::type distance;
                BOOST_MPL_ASSERT((mpl::not_<boost::is_same<iterator, typename mpl::end<PropertyVector>::type > >));

                typedef typename dcm::parser_generator<Property, Sys, iterator_type>::generator generator;

                std::string s;
                std::back_insert_iterator<std::string> out(s);
                generator gen;
                dcm::parser_generator<Property, Sys, iterator_type>::init(gen);
                boost::spirit::karma::generate(out, gen, fusion::at<distance>(sequence));
                stream << s;
            };
        };

        template<typename Derived, typename Sig>
        void generate_properties(std::ostream stream, boost::shared_ptr< dcm::Object<Sys, Derived, Sig> > obj)  {

            typedef typename dcm::Object<Sys, Derived, Sig>::Sequence Properties;
        };

        struct inheriter {

            inheriter() :  m_this((Sys*) this) {}
            Sys* m_this;

            void saveState(std::ostream& stream) {
                typedef std::back_insert_iterator<std::string> iterator_type;
                typedef vertex_container<typename Sys::Cluster> v_container;
                typedef typename Sys::Cluster::global_vertex_iterator viter;

                std::string s;
                std::back_insert_iterator<std::string> out(s);

                karma::rule<iterator_type, boost::iterator_range<viter>(), karma::locals<int> > r1;
                karma::rule<iterator_type, dcm::GlobalVertex()> r2;
                r2 = karma::int_<<'>'<<karma::lit(karma::_val);
                r1 = (karma::lit("<Vertex id=")<<r2<<"</Vertex>") % karma::eol;
                boost::iterator_range<viter> range(m_this->m_cluster.globalVertices().first,m_this->m_cluster.globalVertices().second);

                karma::generate(out, r1, range);
                stream << s;

                std::pair<viter, viter> it = m_this->m_cluster.globalVertices();
                for(; it.first != it.second; it.first++)
                    std::cout<<"this: "<<(int) *(it.first)<<std::endl;
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

#endif //DCM_MODULEPARSER_H




