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

#ifndef DCM_PARSER_H
#define DCM_PARSER_H

#include <iosfwd>

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/support_istream_iterator.hpp>
#include <boost/spirit/include/qi_string.hpp>
#include <boost/spirit/include/phoenix.hpp>

#include "opendcm/core/clustergraph.hpp"

namespace qi = boost::spirit::qi;
namespace ascii = boost::spirit::ascii;
namespace phx = boost::phoenix;

namespace dcm {

typedef boost::spirit::istream_iterator IIterator;

struct sp : qi::grammar<IIterator, std::string()> {
  
  qi::rule<IIterator, std::string()> start;
  sp() : sp::base_type(start) {
    start %= +qi::char_;
  };
  static void print(std::string& s) {
    std::cout<<"parsed string:"<<std::endl<<s<<std::endl<<"done print string"<<std::endl;
  };
};  

template<typename Prop>
struct empty_pop_parser : public qi::grammar<IIterator, typename Prop::type()> {
        qi::rule<IIterator, typename Prop::type()> start;
        empty_pop_parser(): empty_pop_parser::base_type(start) {
	    start = qi::eps(false);
	};
};
    
template<typename Prop, typename Par>
struct prop_parser : qi::grammar<IIterator, typename Prop::type(), qi::space_type> { 
    
    typename Par::parser subrule;
    sp str;
    qi::rule<IIterator, typename Prop::type(), qi::space_type> start;
    prop_parser(): prop_parser<Prop, Par>::base_type(start) {
	Par::init(subrule);
	start =  qi::lit("<Property>") >> subrule >> qi::lit("</Property>");
    };
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

        prop_par(): prop_par<Sys, PropertyList>::base_type(prop) {

	    prop =       -(qi::eps(valid<0>::value) >> fusion::at<index<0> >(rules)[phx::at_c<index<0>::value>(qi::_val) = qi::_1])
		      >> -(qi::eps(valid<1>::value) >> fusion::at<index<1> >(rules)[phx::at_c<index<1>::value>(qi::_val) = qi::_1])
		      >> -(qi::eps(valid<2>::value) >> fusion::at<index<2> >(rules)[phx::at_c<index<2>::value>(qi::_val) = qi::_1])
		      >> -(qi::eps(valid<3>::value) >> fusion::at<index<3> >(rules)[phx::at_c<index<3>::value>(qi::_val) = qi::_1])
		      >> -(qi::eps(valid<4>::value) >> fusion::at<index<4> >(rules)[phx::at_c<index<4>::value>(qi::_val) = qi::_1])
		      >> -(qi::eps(valid<5>::value) >> fusion::at<index<5> >(rules)[phx::at_c<index<5>::value>(qi::_val) = qi::_1])
		      >> -(qi::eps(valid<6>::value) >> fusion::at<index<6> >(rules)[phx::at_c<index<6>::value>(qi::_val) = qi::_1])
		      >> -(qi::eps(valid<7>::value) >> fusion::at<index<7> >(rules)[phx::at_c<index<7>::value>(qi::_val) = qi::_1])
		      >> -(qi::eps(valid<8>::value) >> fusion::at<index<8> >(rules)[phx::at_c<index<8>::value>(qi::_val) = qi::_1])
		      >> -(qi::eps(valid<9>::value) >> fusion::at<index<9> >(rules)[phx::at_c<index<9>::value>(qi::_val) = qi::_1]);
	};
    };
    
    
    
    
    
    
    
    
    
    
    
    
    template<typename Obj, typename Sys>
    struct empty_obj_parser : public qi::grammar<IIterator, boost::shared_ptr<Obj>(Sys*), qi::space_type> {
	    qi::rule<IIterator, boost::shared_ptr<Obj>(Sys*), qi::space_type> start;
	    empty_obj_parser(): empty_obj_parser::base_type(start) {
		//start = qi::eps(false);
	    };
    };
    
        //grammar for a single object
    template<typename Sys, typename Object, typename Par>
    struct obj_parser : public qi::grammar<IIterator, boost::shared_ptr<Object>(Sys*), qi::space_type> {
        typename Par::parser subrule;
        qi::rule<IIterator, boost::shared_ptr<Object>(Sys*), qi::space_type> start;
        prop_par<Sys, typename Object::Sequence > prop;

        obj_parser(): obj_parser::base_type(start) {
	    Par::init(subrule);
	    start = qi::lit("<Object>") >> subrule(qi::_r1)[qi::_val = qi::_1]
		    >> qi::eps(qi::_val)[ phx::bind(&Sys::template push_back<Object>, qi::_r1, qi::_val)]
		    >> prop[phx::bind(&obj_parser::setProperties, qi::_val, qi::_1)]
		    >> qi::lit("</Object>");
	};
        static void setProperties(boost::shared_ptr<Object> ptr, typename details::pts<typename Object::Sequence>::type& seq) {
	  if(ptr) ptr->m_properties = seq;
	};
    };

    //when objects should not be generated we need to get a empy rule, as obj_rule_init
    //trys always to access the rules attribute and when the parser_generator trait is not
    //specialitzed it's impossible to have the attribute type right in the unspecialized trait
    template<typename Sys, typename seq, typename state>
    struct obj_parser_fold : mpl::fold< seq, state,
            mpl::if_< parser_parse<mpl::_2, Sys>,
            mpl::push_back<mpl::_1,
            obj_parser<Sys, mpl::_2, dcm::parser_parser<mpl::_2, Sys, IIterator> > >,
            mpl::push_back<mpl::_1, empty_obj_parser<mpl::_2, Sys> > > > {};

    //currently max. 10 objects are supported
    template<typename Sys>
    struct obj_par : public qi::grammar<IIterator, 
					typename details::sps<typename Sys::objects>::type(Sys*),
					qi::space_type> {

	typedef typename Sys::objects ObjectList;
      
        //create a vector with the appropriate rules for all objects. Do this with the rule init struct, as it gives
        //automatic initialisation of the rules when the objects are created
        typedef typename obj_parser_fold<Sys, ObjectList, mpl::vector<> >::type init_rules_vector;
        //push back a empty rule so that we know where to go when nothing is to do
        typedef typename mpl::push_back<init_rules_vector,
					empty_obj_parser<typename mpl::back<ObjectList>::type, Sys> >::type rules_vector;

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
        qi::rule<IIterator, typename details::sps<ObjectList>::type(Sys*), qi::space_type> obj;

        obj_par(): obj_par<Sys>::base_type(obj) {

	    obj =      -(qi::eps(valid<0>::value) >> fusion::at<index<0> >(rules)(qi::_r1)[phx::at_c<index<0>::value>(qi::_val) = qi::_1])
		    >> -(qi::eps(valid<1>::value) >> fusion::at<index<1> >(rules)(qi::_r1)[phx::at_c<index<1>::value>(qi::_val) = qi::_1])
		    >> -(qi::eps(valid<2>::value) >> fusion::at<index<2> >(rules)(qi::_r1)[phx::at_c<index<2>::value>(qi::_val) = qi::_1])
		    >> -(qi::eps(valid<3>::value) >> fusion::at<index<3> >(rules)(qi::_r1)[phx::at_c<index<3>::value>(qi::_val) = qi::_1])
		    >> -(qi::eps(valid<4>::value) >> fusion::at<index<4> >(rules)(qi::_r1)[phx::at_c<index<4>::value>(qi::_val) = qi::_1])
		    >> -(qi::eps(valid<5>::value) >> fusion::at<index<5> >(rules)(qi::_r1)[phx::at_c<index<5>::value>(qi::_val) = qi::_1])
		    >> -(qi::eps(valid<6>::value) >> fusion::at<index<6> >(rules)(qi::_r1)[phx::at_c<index<6>::value>(qi::_val) = qi::_1])
		    >> -(qi::eps(valid<7>::value) >> fusion::at<index<7> >(rules)(qi::_r1)[phx::at_c<index<7>::value>(qi::_val) = qi::_1])
		    >> -(qi::eps(valid<8>::value) >> fusion::at<index<8> >(rules)(qi::_r1)[phx::at_c<index<8>::value>(qi::_val) = qi::_1])
		    >> -(qi::eps(valid<9>::value) >> fusion::at<index<9> >(rules)(qi::_r1)[phx::at_c<index<9>::value>(qi::_val) = qi::_1]);

	};
    };
    
    
    
    
    
    
    
    
    
    
    
    
    
template<typename Sys>
struct Injector {
  
  void setClusterProperties(typename Sys::Cluster* cluster,
			    typename details::pts<typename Sys::Cluster::cluster_properties>::type& prop) {
    cluster->m_cluster_bundle = prop;
  };
  void setVertexProperties(typename Sys::Cluster* cluster, LocalVertex v,
			    typename details::pts<typename Sys::vertex_properties>::type& prop) {
    fusion::at_c<0>(cluster->operator[](v)) = prop;
  };
  void setVertexObjects(typename Sys::Cluster* cluster, LocalVertex v,
			typename details::sps<typename Sys::objects>::type& obj) {
    fusion::at_c<1>(cluster->operator[](v)) = obj;
  };
  
  void setEdgeProperties(typename Sys::Cluster* cluster, LocalEdge e,
			    typename details::pts<typename Sys::edge_properties>::type& prop) {
    fusion::at_c<0>(cluster->operator[](e)) = prop;
  };
  void addCluster(typename Sys::Cluster* cluster, fusion::vector2<GlobalVertex, typename Sys::Cluster*>& vec) {
      LocalVertex v = cluster->getLocalVertex(fusion::at_c<0>(vec)).first;
      cluster->m_clusters[v] = fusion::at_c<1>(vec);
  };
    
    
};
    
template<typename Sys>
struct parser : qi::grammar<IIterator, Sys*(), qi::space_type> {

    parser(Sys& s);

    qi::rule<IIterator, Sys*(), qi::space_type> start;
    
    qi::rule<IIterator, fusion::vector2<GlobalVertex, typename Sys::Cluster*>(Sys*), qi::space_type> cluster;
    prop_par<Sys, typename Sys::Cluster::cluster_properties> cluster_prop;
    
    qi::rule<IIterator, 
	     fusion::vector2<LocalVertex, GlobalVertex>(typename Sys::Cluster*, Sys*), qi::space_type> vertex;
    prop_par<Sys, typename Sys::vertex_properties> vertex_prop;
    obj_par<Sys> objects;
    
    qi::rule<IIterator, 
	     fusion::vector4<LocalEdge, GlobalEdge, bool, bool>(typename Sys::Cluster*, Sys*),
	     qi::space_type> edge;
    qi::rule<IIterator, typename Sys::Cluster::edge_bundle_single(Sys*), qi::space_type> global_edge;
    prop_par<Sys, typename Sys::edge_properties> edge_prop;

    sp str;
    Sys& system;
    Injector<Sys> in;
};
  

template<typename Sys>
parser<Sys>::parser(Sys& s) : parser<Sys>::base_type(start), system(s) {
  
    vertex = ascii::string("<Vertex")[qi::_val = phx::bind(&Sys::Cluster::addVertex, qi::_r1)] >> ("id=") 
	     >> qi::int_[phx::at_c<1>(qi::_val) = phx::bind(&Sys::Cluster::setGlobalVertex, qi::_r1, phx::at_c<0>(qi::_val), qi::_1)]
	     >> '>' >> vertex_prop[phx::bind(&Injector<Sys>::setVertexProperties, in, qi::_r1, phx::at_c<0>(qi::_val), qi::_1)]
	     >> objects(qi::_r2)[phx::bind(&Injector<Sys>::setVertexObjects, in, qi::_r1, phx::at_c<0>(qi::_val), qi::_1)]
	     >> ("</Vertex>");

    global_edge = qi::lit("<GlobalEdge") >> qi::lit("id=") >> qi::int_[phx::bind(&GlobalEdge::ID, phx::at_c<1>(qi::_val)) = qi::_1] 
		>> qi::lit("source=") >> qi::int_[phx::bind(&GlobalEdge::source, phx::at_c<1>(qi::_val)) = qi::_1]
		>> qi::lit("target=") >> qi::int_[phx::bind(&GlobalEdge::target, phx::at_c<1>(qi::_val)) = qi::_1] >> '>'
		>> objects(qi::_r1)[phx::at_c<0>(qi::_val) = qi::_1] >> "</GlobalEdge>";
		     
    edge   = (qi::lit("<Edge") >> "source=" >> qi::int_ >> "target=" >> qi::int_ >> '>')[qi::_val = phx::bind((&Sys::Cluster::addEdgeGlobal), qi::_r1, qi::_1, qi::_2)] 
	     >> edge_prop[phx::bind(&Injector<Sys>::setEdgeProperties, in, qi::_r1, phx::at_c<0>(qi::_val), qi::_1)]
	     >> *global_edge(qi::_r2)
	     >> ("</Edge>");

    cluster = qi::lit("<Cluster id=")[phx::at_c<1>(qi::_val) = phx::new_<typename Sys::Cluster>()]
              >> qi::int_[phx::at_c<0>(qi::_val) = qi::_1] >> ">" 
              >> cluster_prop[phx::bind(&Injector<Sys>::setClusterProperties, in, phx::at_c<1>(qi::_val), qi::_1)]
	      >> *vertex(phx::at_c<1>(qi::_val), qi::_r1)
	      >> *edge(phx::at_c<1>(qi::_val), qi::_r1)
	      >> *(cluster(qi::_r1)[phx::bind(&Injector<Sys>::addCluster, in, phx::at_c<1>(qi::_val), qi::_1)])
	      >> "</Cluster>";

    start = qi::no_skip[ascii::string("<Cluster id=0>")] >> cluster_prop[phx::bind(&Injector<Sys>::setClusterProperties, in, &phx::bind(&Sys::m_cluster, qi::_val), qi::_1)]
	    >> *(vertex(&phx::bind(&Sys::m_cluster, qi::_val), qi::_val))
	    >> *(edge(&phx::bind(&Sys::m_cluster, qi::_val), qi::_val)) 
	    >> *(cluster(qi::_val)[phx::bind(&Injector<Sys>::addCluster, in, &phx::bind(&Sys::m_cluster, qi::_val), qi::_1)])
	    >> "</Cluster>" >> str[phx::bind(&sp::print, qi::_1)];
};

}
#endif //DCM_PARSER_H
