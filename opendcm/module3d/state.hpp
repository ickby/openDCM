/*
    openDCM, dimensional constraint manager
    Copyright (C) 2013  Stefan Troeger <stefantroeger@gmx.net>

    This library is free software; you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License along
    with this library; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/

#ifndef DCM_MODULE3D_STATE_HPP
#define DCM_MODULE3D_STATE_HPP

#include "module.hpp"
#include <opendcm/moduleState/traits.hpp>
#include <opendcm/core/clustergraph.hpp>

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/karma.hpp>
#include <boost/spirit/include/phoenix.hpp>
#include <boost/phoenix/function/adapt_function.hpp>

#include <ios>

namespace karma = boost::spirit::karma;
namespace qi = boost::spirit::qi;
namespace karma_ascii = boost::spirit::karma::ascii;
namespace phx = boost::phoenix;

namespace dcm {
  
namespace details {
   
template<typename Sys>
struct getModule3D {
    typedef typename system_traits<Sys>::template getModule<m3d>::type type;
};

    
struct geom_visitor : public boost::static_visitor<std::string> {
   
    template<typename T>
    std::string operator()(T& i) const {
      
	//we use stings in case new geometry gets added and the weights shift, meaning: backwards 
	//compatible
        std::string type;
	switch( geometry_traits<T>::tag::weight::value ) {
	  case 0:
	    return "direction";
	  case 1:
	    return "point";
	  case 2:
	    return "line";
	  case 3:
	    return "plane";
	  case 4:
	    return "cylinder";
	  default:
	    return "unknown";
	};
    };
}; 

template<typename T>
std::string getWeight(boost::shared_ptr<T> ptr) {
    geom_visitor v;
    return ptr->apply(v);
};

template<typename T>
struct get_weight {
  typedef typename geometry_traits<T>::tag::weight type;
};

//search the first type in the typevector with the given weight
template<typename Vector, int Weight>
struct getWeightType {
  typedef typename mpl::find_if<Vector, boost::is_same<get_weight<mpl::_1>, mpl::int_<Weight> > >::type iter;
  typedef typename mpl::deref<iter>::type type;
};

template <typename Geom, typename Row, typename Value>
bool VectorOutput(Geom &v, Row& r, Value& val) {
	  
            if (r < v->m_global.rows()) {
	      
                val = v->m_global(r++);
                return true; // output continues
            }
            return false;    // fail the output
};

template <typename Geom, typename Row, typename Value>
bool VectorInput(Geom &v, Row& r, Value& val) {
	    
	    v.conservativeResize(r+1);
            v(r++) = val;
            return true; // output continues
};

template<typename Geom>
struct inject_set {

    template<typename Vec, typename Obj>
    static void apply(Vec& v, Obj& g) {
	Geom gt;
	(typename geometry_traits<Geom>::modell()).template inject<double,
		  typename geometry_traits<Geom>::accessor >(gt, v);
	g->set(gt);
    };
};
//spezialisation if no type in the typelist hs the right weight
template<>
struct inject_set<mpl_::void_> {

    template<typename Obj, typename Vec>
    static void apply(Vec& v, Obj& g) {
      //TODO:throw   
    };
};

template<typename System>
bool Create(System* sys, std::string& type, 
	      boost::shared_ptr<typename details::getModule3D<System>::type::Geometry3D>& geom, 
	      typename System::Kernel::Vector& v) {
  
  typedef typename details::getModule3D<System>::type::geometry_types Typelist;
  
  if(!type.compare("direction") ) {
    inject_set<typename getWeightType<Typelist, 0>::type>::apply(v, geom);
  }  
  else if(!type.compare("point") ) {
    inject_set<typename getWeightType<Typelist, 1>::type>::apply(v, geom);
  }
  else if(!type.compare("line") ) {
    inject_set<typename getWeightType<Typelist, 2>::type>::apply(v, geom);
  }
  else if(!type.compare("plane") ) {
    inject_set<typename getWeightType<Typelist, 3>::type>::apply(v, geom);
  }
  else if(!type.compare("cylinder") ) {
    inject_set<typename getWeightType<Typelist, 4>::type>::apply(v, geom);
  };
  return true;
};

// define a new real number formatting policy
template <typename Num>
struct scientific_policy : karma::real_policies<Num>
{
    // we want the numbers always to be in scientific format
    static int floatfield(Num n) { return std::ios::scientific; }
    static unsigned precision(Num n) {return 16;};
};

// define a new generator type based on the new policy
typedef karma::real_generator<double, scientific_policy<double> > science_type;
static science_type const scientific = science_type();
} //details
} //dcm
 
BOOST_PHOENIX_ADAPT_FUNCTION( bool, vector_out, dcm::details::VectorOutput, 3)
BOOST_PHOENIX_ADAPT_FUNCTION( bool, vector_in,  dcm::details::VectorInput, 3)
BOOST_PHOENIX_ADAPT_FUNCTION( bool, create,  dcm::details::Create, 4)
 
namespace dcm {  
  
template<typename System>
struct parser_generate< typename details::getModule3D<System>::type::Geometry3D , System>
  : public mpl::true_{};

template<typename System, typename iterator>
struct parser_generator< typename details::getModule3D<System>::type::Geometry3D , System, iterator > {

    typedef typename details::getModule3D<System>::type::Geometry3D  Geometry;
    typedef karma::rule<iterator, boost::shared_ptr<Geometry>(), karma::locals<int> > generator;
    static void init(generator& r) {
       r = karma::lit("<type>Geometry3D</type>\n<class>")
	    << karma_ascii::string[karma::_1 = phx::bind(&details::getWeight<Geometry>, karma::_val)]
	    << "</class>" << karma::eol << "<value>" 
	    << (details::scientific[ boost::spirit::_pass = vector_out(karma::_val, karma::_a, karma::_1) ] % ' ')
	    << "</value>";
    };
};


template<typename System>
struct parser_generate< typename details::getModule3D<System>::type::vertex_prop , System>
  : public mpl::true_{};

template<typename System, typename iterator>
struct parser_generator< typename details::getModule3D<System>::type::vertex_prop , System, iterator > {

    typedef karma::rule<iterator, GlobalVertex()> generator;
    static void init(generator& r) {
        r = karma::lit("<type>Vertex</type>")
	    << karma::eol << "<value>" << karma::int_ << "</value>";
    };
};


template<typename System>
struct parser_generate<typename details::getModule3D<System>::type::fix_prop, System> : public mpl::true_ {};

template<typename System, typename iterator>
struct parser_generator<typename details::getModule3D<System>::type::fix_prop, System, iterator> {
    typedef karma::rule<iterator, bool&()> generator;

    static void init(generator& r) {
        r = karma::lit("<type>Fix</type>\n<value>") << karma::bool_ <<"</value>";
    };
};

/****************************************************************************************************/
/****************************************************************************************************/

template<typename System>
struct parser_parse< typename details::getModule3D<System>::type::Geometry3D , System>
  : public mpl::true_{};

template<typename System, typename iterator>
struct parser_parser< typename details::getModule3D<System>::type::Geometry3D, System, iterator > {

    typedef typename details::getModule3D<System>::type::Geometry3D  object_type;
    typedef typename System::Kernel Kernel;
    
    typedef qi::rule<iterator, boost::shared_ptr<object_type>(System*), qi::space_type, qi::locals<std::string, typename Kernel::Vector, int> > parser;
    static void init(parser& r) {
        r = qi::lit("<type>Geometry3D</type>")[ qi::_val =  phx::construct<boost::shared_ptr<object_type> >( phx::new_<object_type>(*qi::_r1))]
		  >> "<class>" >> (+qi::char_("a-zA-Z"))[qi::_a = phx::construct<std::string>(phx::begin(qi::_1), phx::end(qi::_1))] >> "</class>"
		  >> "<value>" >> *qi::double_[ vector_in(qi::_b, qi::_c, qi::_1) ] >> "</value>"
		  >> qi::eps[ create(qi::_r1, qi::_a, qi::_val, qi::_b) ];
    };
};

template<typename System>
struct parser_parse< typename details::getModule3D<System>::type::vertex_prop, System>
  : public mpl::true_{};

template<typename System, typename iterator>
struct parser_parser< typename details::getModule3D<System>::type::vertex_prop, System, iterator > {
  
    typedef qi::rule<iterator, GlobalVertex(), qi::space_type> parser;
    static void init(parser& r) {
        r = qi::lit("<type>Vertex</type>") >> "<value>" >> qi::int_ >> "</value>";
    };
};

template<typename System>
struct parser_parse< typename details::getModule3D<System>::type::fix_prop, System>
  : public mpl::true_{};

template<typename System, typename iterator>
struct parser_parser< typename details::getModule3D<System>::type::fix_prop, System, iterator > {
  
    typedef qi::rule<iterator, bool(), qi::space_type> parser;
    static void init(parser& r) {
        r = qi::lit("<type>Fix</type>") >> "<value>" >> qi::bool_ >> "</value>";
    };
};

}


#endif //DCM_MODULE3D_STATE_HPP