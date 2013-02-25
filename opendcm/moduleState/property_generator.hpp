#ifndef DCM_PROPERTY_GENERATOR_H
#define DCM_PROPERTY_GENERATOR_H

#include <boost/fusion/include/as_vector.hpp>
#include <boost/spirit/include/karma.hpp>

#include "traits.hpp"

using namespace boost::spirit::karma;
namespace fusion = boost::fusion;

namespace dcm {
  
typedef std::ostream_iterator<char> Iterator;

namespace details {
  
      //a grammar that does nothing exept failing
    struct empty_grammar : public grammar<Iterator> {
        rule<Iterator> start;
        empty_grammar(): empty_grammar::base_type(start) {
	    start = eps(false);
	};
    };

    //grammar for a single property
    template<typename Prop, typename Gen>
    struct prop_grammar : public grammar<Iterator, typename Prop::type()> {
        typename Gen::generator subrule;
        rule<Iterator, typename Prop::type()> start;
        prop_grammar();
    };

    template<typename Sys, typename seq, typename state>
    struct prop_generator_fold : mpl::fold< seq, state,
            mpl::if_< parser_generate<mpl::_2, Sys>,
            mpl::push_back<mpl::_1,
            prop_grammar<mpl::_2, dcm::parser_generator<mpl::_2, Sys, Iterator> > >,
            mpl::push_back<mpl::_1, empty_grammar > > > {};

    //grammar for a fusion sequence of properties. currently max. 10 properties are supported
    template<typename Sys, typename PropertyList>
    struct prop_gen : grammar<Iterator, typename details::pts<PropertyList>::type()> {

        //create a vector with the appropriate rules for all properties. Do this with the rule init struct, as it gives
        //automatic initialisation of the rules when the objects are created
        typedef typename prop_generator_fold<Sys, PropertyList,mpl::vector<> >::type init_rules_vector;
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
        rule<Iterator, typename details::pts<PropertyList>::type()> prop;

        prop_gen();
    };
    
    //special prop classes for better externalisaton, therefore the outside constructor to avoid auto inline
    template<typename Sys> 
    struct cluster_prop_gen : public prop_gen<Sys, typename Sys::Cluster::cluster_properties> {
      cluster_prop_gen();
    };
   
    template<typename Sys>     
    struct vertex_prop_gen : public prop_gen<Sys, typename Sys::Cluster::vertex_properties> {
      vertex_prop_gen();
    };
    
    template<typename Sys> 
    struct edge_prop_gen : public prop_gen<Sys, typename Sys::Cluster::edge_properties> {
      edge_prop_gen();
    };

}//details
}//dcm

#ifndef USE_EXTERNAL
  #include "property_generator_imp.hpp"
#endif

#endif //DCM_PROPERTY_GENERATOR_H