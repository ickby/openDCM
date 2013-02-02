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
#include <opendcm/core/clustergraph.hpp>
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
        struct prop_rule_init {
            typename T::generator rule;
            prop_rule_init() {
                T::init(rule);
            };
        };

        //currently max. 10 properties are supported
        template<typename Iterator, typename PropertyList>
        struct prop_gen : karma::grammar<Iterator, typename details::pts<PropertyList>::type()> {

            //create a vector with the appropriate rules for all properties. Do this with the rule init struct, as it gives
            //automatic initialisation of the rules when the objects are created
            typedef typename mpl::transform<PropertyList,
            prop_rule_init< dcm::parser_generator<mpl::_1, Sys, Iterator> > >::type init_rules_vector;
            //add a empty rule to the end so that we can call it everytime our propertvector is smaller 10
            typedef typename mpl::push_back<init_rules_vector, prop_rule_init< dcm::parser_generator<int, Sys, Iterator> > >::type rules_vector;
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
            rule<Iterator, typename details::pts<PropertyList>::type()> prop;

            prop_gen() : prop_gen::base_type(start) {

            prop =   ((karma::eps(valid<0>::value) << fusion::at<index<0> >(rules).rule)[karma::_1 = phx::at_c<index<0>::value >(karma::_val)])
                   | ((karma::eps(valid<1>::value) << (fusion::at<index<1> >(rules).rule)[karma::_1 = phx::at_c<index<1>::value>(karma::_val)])
                   | ((karma::eps(valid<2>::value) << (fusion::at<index<2> >(rules).rule)[karma::_1 = phx::at_c<index<2>::value>(karma::_val)])
                   | ((karma::eps(valid<3>::value) << (fusion::at<index<3> >(rules).rule)[karma::_1 = phx::at_c<index<3>::value>(karma::_val)])
                   | ((karma::eps(valid<4>::value) << (fusion::at<index<4> >(rules).rule)[karma::_1 = phx::at_c<index<4>::value>(karma::_val)])
                   | ((karma::eps(valid<5>::value) << (fusion::at<index<5> >(rules).rule)[karma::_1 = phx::at_c<index<5>::value>(karma::_val)])
                   | ((karma::eps(valid<6>::value) << (fusion::at<index<6> >(rules).rule)[karma::_1 = phx::at_c<index<6>::value>(karma::_val)])
                   | ((karma::eps(valid<7>::value) << (fusion::at<index<7> >(rules).rule)[karma::_1 = phx::at_c<index<7>::value>(karma::_val)])
                   | ((karma::eps(valid<8>::value) << (fusion::at<index<8> >(rules).rule)[karma::_1 = phx::at_c<index<8>::value>(karma::_val)])
                   | (karma::eps(valid<9>::value) << (fusion::at<index<9> >(rules).rule)[karma::_1 = phx::at_c<index<9>::value>(karma::_val)])))))))));

            start = &prop[karma::_1=karma::_val] << karma::lit("<Property>") << karma::eol << prop[karma::_1=karma::_val] << karma::eol << karma::lit("</Property>");
        };
        };

        template<typename T>
        struct getSequence {
            typedef typename T::Object::Sequence type;
        };

        template<typename Iterator, typename Object, typename Gen>
        struct obj_rule_init {
            typename Gen::generator subrule,rule;
            obj_rule_init() {
                Gen::init(subrule);

                prop_gen<Iterator, typename Object::Sequence > prop;
                rule = subrule << prop[phx::bind(&obj_rule_init::getProperties, this, karma::_val, karma::_1)];
            };
            void getProperties(boost::shared_ptr<Object> ptr, typename details::pts<typename Object::Sequence>::type& seq) {
                seq = ptr->m_properties;
            };
        };
        template<typename Iterator>
        struct empty_rule_init {
            karma::rule<Iterator> rule;
        };

        //when objects should not be generated we need to get a empy rule, as obj_rule_init
        //trys always to access the rules attribute and when the parser_generator trait is not
        //specialitzed it's impossible to have the attribute type right in the unspecialized trait
        template<typename seq, typename state, typename Iterator>
        struct obj_generator_fold : mpl::fold< seq, state,
                mpl::if_< parser_generate<mpl::_2, Sys>,
                mpl::push_back<mpl::_1,
                obj_rule_init<Iterator, mpl::_2, dcm::parser_generator<mpl::_2, Sys, Iterator> > >,
                mpl::push_back<mpl::_1, empty_rule_init<Iterator> > > > {};

        //currently max. 10 properties are supported
        template<typename Iterator, typename ObjectList>
        struct obj_gen : karma::grammar<Iterator, typename details::sps<ObjectList>::type()> {

            //create a vector with the appropriate rules for all objects. Do this with the rule init struct, as it gives
            //automatic initialisation of the rules when the objects are created
            typedef typename obj_generator_fold<ObjectList,mpl::vector<>, Iterator >::type init_rules_vector;
            //push back a empty rule so that we know where to go when nothing is to do
            typedef typename mpl::push_back<init_rules_vector, empty_rule_init<Iterator> >::type rules_vector;

            //create the fusion sequence of our rules
            typedef typename fusion::result_of::as_vector<rules_vector>::type rules_sequnce;

            //this struct returns the right accessvalue for the sequences. If we access a value bigger than the property vector size
            //we use the last rule, as we made sure this is an empty one
            template<int I>
        struct index : public mpl::if_< mpl::less<mpl::int_<I>, mpl::size<ObjectList> >,
            mpl::int_<I>, typename mpl::size<ObjectList>::prior >::type {};
            //this struct tells us if we should execute the generator
            template<int I>
        struct valid : public mpl::less< mpl::int_<I>, mpl::size<ObjectList> > {};

            rules_sequnce rules;
            rule<Iterator, typename details::sps<ObjectList>::type()> start;
            rule<Iterator, typename details::sps<ObjectList>::type()> obj;

            obj_gen() : obj_gen::base_type(start) {

            obj =    (karma::eps(valid<0>::value) << karma::eps(phx::at_c<index<0>::value>(karma::_val)) << fusion::at<index<0> >(rules).rule[karma::_1 = phx::at_c<index<0>::value>(karma::_val)])
                  | ((karma::eps(valid<1>::value) << karma::eps(phx::at_c<index<1>::value>(karma::_val)) << (fusion::at<index<1> >(rules).rule)[karma::_1 = phx::at_c<index<1>::value>(karma::_val)])
                  | ((karma::eps(valid<2>::value) << karma::eps(phx::at_c<index<2>::value>(karma::_val)) << (fusion::at<index<2> >(rules).rule)[karma::_1 = phx::at_c<index<2>::value>(karma::_val)])
                  | ((karma::eps(valid<3>::value) << karma::eps(phx::at_c<index<3>::value>(karma::_val)) << (fusion::at<index<3> >(rules).rule)[karma::_1 = phx::at_c<index<3>::value>(karma::_val)])
                  | ((karma::eps(valid<4>::value) << karma::eps(phx::at_c<index<4>::value>(karma::_val)) << (fusion::at<index<4> >(rules).rule)[karma::_1 = phx::at_c<index<4>::value>(karma::_val)])
                  | ((karma::eps(valid<5>::value) << karma::eps(phx::at_c<index<5>::value>(karma::_val)) << (fusion::at<index<5> >(rules).rule)[karma::_1 = phx::at_c<index<5>::value>(karma::_val)])
                  | ((karma::eps(valid<6>::value) << karma::eps(phx::at_c<index<6>::value>(karma::_val)) << (fusion::at<index<6> >(rules).rule)[karma::_1 = phx::at_c<index<6>::value>(karma::_val)])
                  | ((karma::eps(valid<7>::value) << karma::eps(phx::at_c<index<7>::value>(karma::_val)) << (fusion::at<index<7> >(rules).rule)[karma::_1 = phx::at_c<index<7>::value>(karma::_val)])
                  | ((karma::eps(valid<8>::value) << karma::eps(phx::at_c<index<8>::value>(karma::_val)) << (fusion::at<index<8> >(rules).rule)[karma::_1 = phx::at_c<index<8>::value>(karma::_val)])
                  | (karma::eps(valid<9>::value) << karma::eps(phx::at_c<index<9>::value>(karma::_val)) << (fusion::at<index<9> >(rules).rule)[karma::_1 = phx::at_c<index<9>::value>(karma::_val)])))))))));

            start = &obj[karma::_1=karma::_val] << (karma::lit("<Object>") << karma::eol << obj[karma::_1=karma::_val] << karma::eol << karma::lit("</Object>"));
        };
        };

        //we need this layer of indirection as we cant access system typedefs in the function declaration of th einheriter
        struct SequenceExtractor {
            Sys* sys;
            SequenceExtractor(Sys* s) : sys(s) {};
            void getVertexProperties(LocalVertex v, typename details::pts< typename Sys::vertex_properties >::type& seq) {
                seq = fusion::at_c<0>(sys->m_cluster[v]);
            };
            void getVertexObjects(LocalVertex v, typename details::sps< typename Sys::objects >::type& seq) {
                seq = fusion::at_c<1>(sys->m_cluster[v]);
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
                obj_gen<iterator_type, typename Sys::objects> obj_gen;
                SequenceExtractor ex(m_this);

                id = karma::int_[phx::bind(&Sys::getVertexID, this, karma::_val, karma::_1)];
                vertex = id << ">\n"
                         << vertex_prop[phx::bind(&SequenceExtractor::getVertexProperties, ex, karma::_val, karma::_1)]
                         << karma::eol
                         << obj_gen[phx::bind(&SequenceExtractor::getVertexObjects, ex, karma::_val, karma::_1)];

                r1 = ((karma::lit("<Vertex id=")<<vertex<<karma::eol<<"</Vertex>")) % karma::eol;

                std::pair<viter,viter> res = boost::vertices(m_this->m_cluster);
                boost::iterator_range<viter> range(res.first,res.second);

                karma::generate(out, r1, range);
                stream << s;

            };

            void loadState(std::istream& stream) {

            };

        private:
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





