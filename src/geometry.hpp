/*
    openDCM, dimensional constraint manager
    Copyright (C) 2012  Stefan Troeger <stefantroeger@gmx.net>

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


#ifndef NS2_GEOMETRY_H
#define NS2_GEOMETRY_H

#include <iostream>

#include <Eigen/Core>

#include <boost/type_traits.hpp>
#include <boost/mpl/assert.hpp>
#include "clustergraph.hpp"
#include "variant.hpp"

namespace dcm {

  struct undefined{};
  struct tag_point{};
  struct tag_line{};
  struct tag_plane{};

  template< typename T>
   struct geometry_traits {
     typedef undefined tag;
   };

namespace detail {
  
  template< typename T>
  struct storage {
    typedef undefined type;
  };
  
  template<>
  struct storage<tag_point> {
    typedef Eigen::Vector3d type;
  };
}
  
template< typename T >
class Geometry {
  
  typedef typename geometry_traits<T>::tag 	tag;
  typedef typename detail::storage<tag>::type 	Storage;
  
  BOOST_MPL_ASSERT(( boost::mpl::not_< boost::is_same<tag, undefined> > ));
  BOOST_MPL_ASSERT(( boost::mpl::not_< boost::is_same<Storage, undefined> >));

public:
    Geometry(T geom) {set(geom);};
    
    T& 		get() {return m_geom;};
    void	set(T geom) {m_geom = geom;};
    
//private:
  
   Storage 		m_storage;
   T 			m_geom;
   Vertex 		m_vertex;
   ClusterGraph*	m_cluster;
   
};

template< typename Tlist> 
class GVariant : public Variant<Tlist> {
  
public:  
    GVariant() : Variant<Tlist>(), Variant<Tlist>::m_uiSelect(0) {};
  
    template< typename T>    
    GVariant(boost::shared_ptr< Geometry<T> > geom) {
       Variant<Tlist>::template assign< boost::shared_ptr< Geometry<T> >, T >(geom);
    };
  
    template<typename T>
    typename T::result_type apply(T& visitor) {
      Variant<Tlist>::template apply_visitor<boost::mpl::int_<0>, ApplyCallback, T>(visitor, boost::mpl::true_());
    };
    
private:
  
    struct ApplyCallback {
     
      template<typename IndexType, typename VisitorType>
      typename VisitorType::result_type operator()(VisitorType& rVisitor, boost::mpl::true_){

	return typename VisitorType::result_type();
      }
      
      template< typename IndexType, typename VisitorType>
      typename VisitorType::result_type operator()(VisitorType& rVisitor, boost::mpl::false_){
	BOOST_ASSERT(false);
	return typename VisitorType::result_type();
      }
    };
};

}

#endif // NS2_GEOMETRY_H
