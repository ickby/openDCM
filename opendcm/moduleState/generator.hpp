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

#ifndef DCM_GENERATOR_H
#define DCM_GENERATOR_H

#include <iosfwd>
#include <boost/shared_ptr.hpp>

#include <opendcm/core/object.hpp>
#include <opendcm/core/property.hpp>
#include <opendcm/core/clustergraph.hpp>
#include <opendcm/modulePart/module.hpp>

#include "traits.hpp"
#include "indent.hpp"
#include <stdio.h>
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
#include <boost/spirit/home/support/container.hpp>

using namespace boost::spirit::karma;

namespace karma = boost::spirit::karma;
namespace phx = boost::phoenix;
namespace fusion = boost::fusion;

namespace dcm {

template<typename Iterator, typename Sys>
struct generator : grammar<Iterator, typename Sys::Cluster&()> {

    typedef typename boost::graph_traits<typename Sys::Cluster>::vertex_iterator viter;
    typedef typename boost::graph_traits<typename Sys::Cluster>::edge_iterator eiter;

    //we need this layer of indirection as we cant access system typedefs in the function declaration of th einheriter
    struct Extractor {
        Sys& sys;
        Extractor(Sys& s) : sys(s) {};
        void getVertexProperties(LocalVertex v, typename details::pts< typename Sys::vertex_properties >::type& seq) {
            seq = fusion::at_c<0>(sys.m_cluster[v]);
        };
        void getEdgeProperties(LocalEdge e, typename details::pts< typename Sys::edge_properties >::type& seq) {
            seq = fusion::at_c<0>(sys.m_cluster[e]);
        };
        void getVertexObjects(LocalVertex v, typename details::sps< typename Sys::objects >::type& seq) {
            seq = fusion::at_c<1>(sys.m_cluster[v]);
        };
	void getEdgeObjects(GlobalEdge v, typename details::sps< typename Sys::objects >::type& seq) {
            seq = fusion::at_c<1>(sys.m_cluster[v]);
        };
        void getVertexID(LocalVertex v, int& id) {
            id = sys.m_cluster.getGlobalVertex(v);
        };
        void getVertexRange(typename Sys::Cluster& cluster, std::vector<LocalVertex>& range) {
            std::pair<viter, viter> res = boost::vertices(cluster);
	    for(;res.first != res.second; res.first++)
	      range.push_back(*res.first);
        };
        void getEdgeRange(typename Sys::Cluster& cluster, std::vector<LocalEdge>& range) {
            //std::pair<eiter, eiter> res = boost::edges(cluster);
	    std::pair<eiter, eiter> res = boost::edges(sys.m_cluster);
	    for(;res.first != res.second; res.first++)
	      range.push_back(*res.first);
        };
        void getEdgeSource(LocalEdge e, int& source) {
	    LocalVertex v = boost::source(e,sys.m_cluster);
            getVertexID(e.m_source, source);
        };
        void getEdgeTarget(LocalEdge e, int& target) {
            getVertexID(boost::target(e,sys.m_cluster), target);
        };
    };

//a grammar that does nothing
    struct empty_grammar : public grammar<Iterator> {
        rule<Iterator> start;
        empty_grammar() : empty_grammar::base_type(start) {};
    };

//grammar for a single property
    template<typename Prop, typename Gen>
    struct prop_grammar : public grammar<Iterator, typename Prop::type()> {
        typename Gen::generator subrule;
        rule<Iterator, typename Prop::type()> start;
        prop_grammar() : prop_grammar::base_type(start) {
            Gen::init(subrule);
            start =  lit("<Property>") << '+' << eol << subrule
                     << '-' << eol << lit("</Property>");
        };
    };

//grammar for a fusion sequence of properties. currently max. 10 properties are supported
    template<typename PropertyList>
    struct prop_gen : grammar<Iterator, typename details::pts<PropertyList>::type()> {

        //create a vector with the appropriate rules for all properties. Do this with the rule init struct, as it gives
        //automatic initialisation of the rules when the objects are created
        typedef typename mpl::transform<PropertyList,
                prop_grammar<mpl::_1, dcm::parser_generator<mpl::_1, Sys, Iterator> > >::type init_rules_vector;
        //add a empty rule to the end so that we can call it everytime our propertvector is smaller 10
        typedef typename mpl::push_back<init_rules_vector, empty_grammar >::type rules_vector;
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

            prop = (eps(valid<0>::value) << fusion::at<index<0> >(rules)[karma::_1 = phx::at_c<index<0>::value >(_val)])
                   | ((eps(valid<1>::value) << (fusion::at<index<1> >(rules))[karma::_1 = phx::at_c<index<1>::value>(_val)])
                      | ((eps(valid<2>::value) << (fusion::at<index<2> >(rules))[karma::_1 = phx::at_c<index<2>::value>(_val)])
                         | ((eps(valid<3>::value) << (fusion::at<index<3> >(rules))[karma::_1 = phx::at_c<index<3>::value>(_val)])
                            | ((eps(valid<4>::value) << (fusion::at<index<4> >(rules))[karma::_1 = phx::at_c<index<4>::value>(_val)])
                               | ((eps(valid<5>::value) << (fusion::at<index<5> >(rules))[karma::_1 = phx::at_c<index<5>::value>(_val)])
                                  | ((eps(valid<6>::value) << (fusion::at<index<6> >(rules))[karma::_1 = phx::at_c<index<6>::value>(_val)])
                                     | ((eps(valid<7>::value) << (fusion::at<index<7> >(rules))[karma::_1 = phx::at_c<index<7>::value>(_val)])
                                        | ((eps(valid<8>::value) << (fusion::at<index<8> >(rules))[karma::_1 = phx::at_c<index<8>::value>(_val)])
                                           | (eps(valid<9>::value) << (fusion::at<index<9> >(rules))[karma::_1 = phx::at_c<index<9>::value>(_val)])))))))));

            start = &prop[karma::_1=_val] << prop[karma::_1=_val];
        };
    };

    template<typename T>
    struct getSequence {
        typedef typename T::Object::Sequence type;
    };

    template<typename Object, typename Gen>
    struct obj_grammar : public grammar<Iterator, boost::shared_ptr<Object>()> {
        typename Gen::generator subrule;
        rule<Iterator, boost::shared_ptr<Object>()> start;
        prop_gen<typename Object::Sequence > prop;

        obj_grammar() : obj_grammar::base_type(start) {
            Gen::init(subrule);
            start = lit("<Object>") << '+' << eol << subrule
                    << -(&prop[phx::bind(&obj_grammar::getProperties, _val, karma::_1)] <<
                         eol << prop[phx::bind(&obj_grammar::getProperties, _val, karma::_1)])
                    << '-' << eol << lit("</Object>");
        };
        static void getProperties(boost::shared_ptr<Object> ptr, typename details::pts<typename Object::Sequence>::type& seq) {
            if(ptr) seq = ptr->m_properties;
            else {
                //TODO: throw
            };
        };
    };

//when objects should not be generated we need to get a empy rule, as obj_rule_init
//trys always to access the rules attribute and when the parser_generator trait is not
//specialitzed it's impossible to have the attribute type right in the unspecialized trait
    template<typename seq, typename state>
    struct obj_generator_fold : mpl::fold< seq, state,
            mpl::if_< parser_generate<mpl::_2, Sys>,
            mpl::push_back<mpl::_1,
            obj_grammar<mpl::_2, dcm::parser_generator<mpl::_2, Sys, Iterator> > >,
            mpl::push_back<mpl::_1, empty_grammar > > > {};

//currently max. 10 properties are supported
    template<typename ObjectList>
    struct obj_gen : public grammar<Iterator, typename details::sps<ObjectList>::type()> {

        //create a vector with the appropriate rules for all objects. Do this with the rule init struct, as it gives
        //automatic initialisation of the rules when the objects are created
        typedef typename obj_generator_fold<ObjectList,mpl::vector<> >::type init_rules_vector;
        //push back a empty rule so that we know where to go when nothing is to do
        typedef typename mpl::push_back<init_rules_vector, empty_grammar >::type rules_vector;

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

            obj = (eps(valid<0>::value) << eps(phx::at_c<index<0>::value>(_val)) << fusion::at<index<0> >(rules)[karma::_1 = phx::at_c<index<0>::value>(_val)])
                  | ((eps(valid<1>::value) << eps(phx::at_c<index<1>::value>(_val)) << (fusion::at<index<1> >(rules))[karma::_1 = phx::at_c<index<1>::value>(_val)])
                     | ((eps(valid<2>::value) << eps(phx::at_c<index<2>::value>(_val)) << (fusion::at<index<2> >(rules))[karma::_1 = phx::at_c<index<2>::value>(_val)])
                        | ((eps(valid<3>::value) << eps(phx::at_c<index<3>::value>(_val)) << (fusion::at<index<3> >(rules))[karma::_1 = phx::at_c<index<3>::value>(_val)])
                           | ((eps(valid<4>::value) << eps(phx::at_c<index<4>::value>(_val)) << (fusion::at<index<4> >(rules))[karma::_1 = phx::at_c<index<4>::value>(_val)])
                              | ((eps(valid<5>::value) << eps(phx::at_c<index<5>::value>(_val)) << (fusion::at<index<5> >(rules))[karma::_1 = phx::at_c<index<5>::value>(_val)])
                                 | ((eps(valid<6>::value) << eps(phx::at_c<index<6>::value>(_val)) << (fusion::at<index<6> >(rules))[karma::_1 = phx::at_c<index<6>::value>(_val)])
                                    | ((eps(valid<7>::value) << eps(phx::at_c<index<7>::value>(_val)) << (fusion::at<index<7> >(rules))[karma::_1 = phx::at_c<index<7>::value>(_val)])
                                       | ((eps(valid<8>::value) << eps(phx::at_c<index<8>::value>(_val)) << (fusion::at<index<8> >(rules))[karma::_1 = phx::at_c<index<8>::value>(_val)])
                                          | (eps(valid<9>::value) << eps(phx::at_c<index<9>::value>(_val)) << (fusion::at<index<9> >(rules))[karma::_1 = phx::at_c<index<9>::value>(_val)])))))))));

            start = &obj[karma::_1=_val] << obj[karma::_1=_val];
        };
    };

    generator(Sys& s) : generator::base_type(start), system(s), ex(s) {

        vertex = int_[phx::bind(&Extractor::getVertexID, ex, _val, karma::_1)] << ">+\n"
                 << -(vertex_prop[phx::bind(&Extractor::getVertexProperties, ex, _val, karma::_1)])
                 << -buffer[ eol << objects[phx::bind(&Extractor::getVertexObjects, ex, _val, karma::_1)] ];

        vertex_range = *(lit("<Vertex id=") << vertex << '-' << eol << lit("</Vertex>") << eol);

        edge =  lit("source=")<<int_[phx::bind(&Extractor::getEdgeSource, ex, _val, karma::_1)]
                <<" target="<<int_[phx::bind(&Extractor::getEdgeTarget, ex, _val, karma::_1)] << ">+\n"
                << -(edge_prop[phx::bind(&Extractor::getEdgeProperties, ex, _val, karma::_1)]);
		//<< -buffer[ eol << objects[phx::bind(&Extractor::getEdgeObjects, ex, _val, karma::_1)] ];

        edge_range = *(lit("<Edge ") << edge << '-' << eol << lit("</Edge>") << eol);

        start = vertex_range[phx::bind(&Extractor::getVertexRange, ex, _val, karma::_1)]
                << edge_range[phx::bind(&Extractor::getEdgeRange, ex, _val, karma::_1)];
    };

    static void test(LocalEdge e1, LocalEdge& e2) {
      std::cout<<"test: "<<e1.m_source<<std::endl;
    };
    rule<Iterator, LocalEdge()> testrule;

    rule<Iterator, typename Sys::Cluster& ()> start;

    rule<Iterator, std::vector<LocalVertex>()> vertex_range;
    rule<Iterator, LocalVertex()> vertex;
    prop_gen<typename Sys::vertex_properties> vertex_prop;

    rule<Iterator, std::vector<LocalEdge>()> edge_range;
    rule<Iterator, LocalEdge()> edge;
    prop_gen<typename Sys::edge_properties> edge_prop;

    obj_gen<typename Sys::objects> objects;
    Extractor ex;

    Sys& system;
};

}//namespace dcm

#endif //DCM_GENERATOR_H

