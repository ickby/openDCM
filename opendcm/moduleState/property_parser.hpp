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

#ifndef DCM_PROPERTY_PARSER_H
#define DCM_PROPERTY_PARSER_H

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/support_istream_iterator.hpp>
#include <boost/spirit/include/qi_string.hpp>
#include <boost/spirit/include/phoenix.hpp>

#include <boost/mpl/less.hpp>
#include <boost/mpl/int.hpp>

namespace fusion = boost::fusion;
namespace qi = boost::spirit::qi;
namespace ascii = boost::spirit::ascii;
namespace phx = boost::phoenix;
namespace mpl = boost::mpl;

namespace dcm {

typedef boost::spirit::istream_iterator IIterator;

namespace details {

template<typename Prop>
struct empty_pop_parser : public qi::grammar<IIterator, typename Prop::type()> {
    qi::rule<IIterator, typename Prop::type()> start;
    empty_pop_parser();
};

template<typename Prop, typename Par>
struct prop_parser : qi::grammar<IIterator, typename Prop::type(), qi::space_type> {

    typename Par::parser subrule;
    qi::rule<IIterator, typename Prop::type(), qi::space_type> start;
    prop_parser();
};

template<typename Sys, typename seq, typename state>
struct prop_parser_fold : mpl::fold< seq, state,
        mpl::if_< dcm::parser_parse<mpl::_2, Sys>,
        mpl::push_back<mpl::_1,
        prop_parser<mpl::_2, dcm::parser_parser<mpl::_2, Sys, IIterator> > >,
        mpl::push_back<mpl::_1, empty_pop_parser<mpl::_2> > > > {};

//grammar for a fusion sequence of properties. currently max. 10 properties are supported
template<typename Sys, typename PropertyList>
struct prop_par : qi::grammar<IIterator, typename details::pts<PropertyList>::type(), qi::space_type> {

    //create a vector with the appropriate rules for all properties. Do this with the rule init struct, as it gives
    //automatic initialisation of the rules when the objects are created
    typedef typename prop_parser_fold<Sys, PropertyList, mpl::vector<> >::type init_rules_vector;
    //add a empty rule to the end so that we can call it everytime our propertvector is smaller 10
    typedef typename mpl::push_back<init_rules_vector, empty_pop_parser<typename mpl::back<PropertyList>::type> >::type rules_vector;
    //create the fusion sequence of our rules
    typedef typename fusion::result_of::as_vector<rules_vector>::type rules_sequnce;

    //this struct returns the right accessvalue for the sequences. If we access a value bigger than the property vector size
    //we use the last rule, as we made sure this is an empty one
    template<int I>
    struct index : public mpl::if_< mpl::less<mpl::int_<I>, mpl::size<PropertyList> >,
    mpl::int_<I>, typename mpl::prior<mpl::size<PropertyList> >::type >::type {};
    //this struct tells us if we should execute the generator
    template<int I>
    struct valid : public mpl::less< mpl::int_<I>, mpl::size<PropertyList> > {};

    rules_sequnce rules;
    qi::rule<IIterator, typename details::pts<PropertyList>::type(), qi::space_type> prop;

    prop_par();
};

    //special prop classes for better externalisaton, therefore the outside constructor to avoid auto inline
    template<typename Sys> 
    struct cluster_prop_par : public prop_par<Sys, typename Sys::Cluster::cluster_properties> {
      cluster_prop_par();
    };
   
    template<typename Sys>     
    struct vertex_prop_par : public prop_par<Sys, typename Sys::Cluster::vertex_properties> {
      vertex_prop_par();
    };
    
    template<typename Sys> 
    struct edge_prop_par : public prop_par<Sys, typename Sys::Cluster::edge_properties> {
      edge_prop_par();
    };

} //DCM
} //details

#ifndef USE_EXTERNAL
  #include "property_parser_imp.hpp"
#endif

#endif
