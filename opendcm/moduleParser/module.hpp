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
#include <opendcm/modulePart/module.hpp>
#include "traits.hpp"
#include "container.hpp"

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

struct ModuleParser {

    template<typename Sys>
    struct type {

        //structure to store a rule and initialize it on construction. Temlate parameter needs to be of type parser_generator.
        template<typename T>
        struct rule_init {
            typename T::generator rule;
            rule_init() {
                T::init(rule);
            };
        };

        //currently max. 10 properties are supported
        template<typename Iterator, typename PropertyList>
        struct prop_gen : karma::grammar<Iterator, typename details::pts<PropertyList>::type()> {

            //create a vector with the appropriate rules for all properties. Do this with the rule init struct, as it gives
            //automatic initialisation of the rules when the objects are created
            typedef typename mpl::transform<PropertyList,
            rule_init< dcm::parser_generator<mpl::_1, Sys, Iterator> > >::type init_rules_vector;
            //add a empty rule to the end so that we can call it everytime our propertvector is smaller 10
            typedef typename mpl::push_back<init_rules_vector, rule_init< dcm::parser_generator<int, Sys, Iterator> > >::type rules_vector;
            //create the fusion sequence of our rules
            typedef typename fusion::result_of::as_vector<rules_vector>::type rules_sequnce;

            //this struct returns the right accessvalue for the sequences. If we access a value bigger than the property vector size
            //we use the last rule, as we made sure this is an empty one
            template<int I>
        struct index : public mpl::if_< mpl::less<mpl::int_<I>, mpl::size<PropertyList> >,
            mpl::int_<I>, mpl::size<PropertyList> >::type {};
            //this struct tells us if we should execute the generator
            template<int I>
        struct valid : public mpl::less< mpl::int_<I>, mpl::size<PropertyList> > {};

            rules_sequnce rules;
            rule<Iterator, typename details::pts<PropertyList>::type()> start;

            prop_gen() : prop_gen::base_type(start) {

            start = -(fusion::at<index<0> >(rules).rule)[karma::_1 = phx::at_c<index<0>::value >(karma::_val)]
                    << -(karma::eps(valid<1>::value) << (fusion::at<index<1> >(rules).rule)[karma::_1 = phx::at_c<index<1>::value >(karma::_val)])
                    << -(karma::eps(valid<2>::value) << (fusion::at<index<2> >(rules).rule)[karma::_1 = phx::at_c<index<2>::value >(karma::_val)])
                    << -(karma::eps(valid<3>::value) << (fusion::at<index<3> >(rules).rule)[karma::_1 = phx::at_c<index<3>::value >(karma::_val)])
                    << -(karma::eps(valid<4>::value) << (fusion::at<index<4> >(rules).rule)[karma::_1 = phx::at_c<index<4>::value >(karma::_val)])
                    << -(karma::eps(valid<5>::value) << (fusion::at<index<5> >(rules).rule)[karma::_1 = phx::at_c<index<5>::value >(karma::_val)])
                    << -(karma::eps(valid<6>::value) << (fusion::at<index<6> >(rules).rule)[karma::_1 = phx::at_c<index<6>::value >(karma::_val)])
                    << -(karma::eps(valid<7>::value) << (fusion::at<index<7> >(rules).rule)[karma::_1 = phx::at_c<index<7>::value >(karma::_val)])
                    << -(karma::eps(valid<8>::value) << (fusion::at<index<8> >(rules).rule)[karma::_1 = phx::at_c<index<8>::value >(karma::_val)])
                    << -(karma::eps(valid<9>::value) << (fusion::at<index<9> >(rules).rule)[karma::_1 = phx::at_c<index<9>::value >(karma::_val)]);
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
                typedef typename boost::graph_traits<typename Sys::Cluster>::vertex_iterator viter;

                std::string s;
                std::back_insert_iterator<std::string> out(s);

                karma::rule<iterator_type, boost::iterator_range<viter>(), karma::locals<LocalVertex> > r1;
                karma::rule<iterator_type, LocalVertex()> id;
                karma::rule<iterator_type, LocalVertex()> vertex;
                prop_gen<iterator_type, typename Sys::vertex_properties> vertex_prop;
                prop_gen<iterator_type, typename Sys::edge_properties> edge_prop;
                SequenceExtractor ex(m_this);

                id = karma::int_[phx::bind(&Sys::getVertexID, this, karma::_val, karma::_1)];
                vertex = id << ">\n" << vertex_prop[phx::bind(&SequenceExtractor::getVertexProperties, ex, karma::_val, karma::_1)] ;

                r1 = ((karma::lit("<Vertex id=")<<vertex<<"</Vertex>")) % karma::eol;

                std::pair<viter,viter> res = boost::vertices(m_this->m_cluster);
                boost::iterator_range<viter> range(res.first,res.second);

                karma::generate(out, r1, range);
                stream << s;

            };

            void loadState(std::istream& stream) {

            };

        private:
            //we need this layer of indirection as we cant access system typedefs in the function declaration of th einheriter
            struct SequenceExtractor {
                Sys* sys;
                SequenceExtractor(Sys* s) : sys(s) {};
                void getVertexProperties(LocalVertex v, typename details::pts< typename Sys::vertex_properties >::type& seq) {
                    seq = fusion::at_c<0>(sys->m_cluster[v]);
                };
            };
            void getVertexID(LocalVertex v, int& id) {
                id = m_this->m_cluster.getGlobalVertex(v);
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




