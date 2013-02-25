#ifndef DCM_PROPERTY_GENERATOR_IMP_H
#define DCM_PROPERTY_GENERATOR_IMP_H

#include "property_generator.hpp"
#include "traits_impl.hpp"
#include <boost/spirit/include/phoenix.hpp>

namespace phx = boost::phoenix;

namespace dcm {

typedef std::ostream_iterator<char> Iterator;

namespace details {

//grammar for a single property
template<typename Prop, typename Gen>
prop_grammar<Prop, Gen>::prop_grammar() : prop_grammar<Prop, Gen>::base_type(start) {
    Gen::init(subrule);
    start =  lit("\n<Property>") << '+' << eol << subrule
             << '-' << eol << lit("</Property>");
};

template<typename Sys, typename PropertyList>
prop_gen<Sys, PropertyList>::prop_gen() : prop_gen<Sys, PropertyList>::base_type(prop) {

    prop =    -(eps(valid<0>::value) << fusion::at<index<0> >(rules)[karma::_1 = phx::at_c<index<0>::value >(_val)])
              << -(eps(valid<1>::value) << fusion::at<index<1> >(rules)[karma::_1 = phx::at_c<index<1>::value>(_val)])
              << -(eps(valid<2>::value) << fusion::at<index<2> >(rules)[karma::_1 = phx::at_c<index<2>::value>(_val)])
              << -(eps(valid<3>::value) << fusion::at<index<3> >(rules)[karma::_1 = phx::at_c<index<3>::value>(_val)])
              << -(eps(valid<4>::value) << fusion::at<index<4> >(rules)[karma::_1 = phx::at_c<index<4>::value>(_val)])
              << -(eps(valid<5>::value) << fusion::at<index<5> >(rules)[karma::_1 = phx::at_c<index<5>::value>(_val)])
              << -(eps(valid<6>::value) << fusion::at<index<6> >(rules)[karma::_1 = phx::at_c<index<6>::value>(_val)])
              << -(eps(valid<7>::value) << fusion::at<index<7> >(rules)[karma::_1 = phx::at_c<index<7>::value>(_val)])
              << -(eps(valid<8>::value) << fusion::at<index<8> >(rules)[karma::_1 = phx::at_c<index<8>::value>(_val)])
              << -(eps(valid<9>::value) << fusion::at<index<9> >(rules)[karma::_1 = phx::at_c<index<9>::value>(_val)]);
};

template<typename Sys>
cluster_prop_gen<Sys>::cluster_prop_gen() : prop_gen<Sys, typename Sys::Cluster::cluster_properties>() {};

template<typename Sys>
vertex_prop_gen<Sys>::vertex_prop_gen() : prop_gen<Sys, typename Sys::Cluster::vertex_properties>() {};

template<typename Sys>
edge_prop_gen<Sys>::edge_prop_gen() : prop_gen<Sys, typename Sys::Cluster::edge_properties>() {};

}//details
}//dcm

#endif //DCM_PROPERTY_GENERATOR_H
