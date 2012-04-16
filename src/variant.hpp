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

#ifndef NS2_VARIANT_H
#define NS2_VARIANT_H

#include <boost/type_traits/is_same.hpp>
#include <boost/type_traits/is_base_of.hpp>
#include <boost/mpl/is_sequence.hpp>
#include <boost/mpl/less.hpp>
#include <boost/mpl/size.hpp>
#include <boost/mpl/at.hpp>
#include <boost/mpl/find.hpp>
#include <boost/mpl/assert.hpp>
#include <boost/static_assert.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/preprocessor/iteration/local.hpp>

#include "geometry.hpp"

namespace dcm  {

template< typename Tlist >
class Variant  {
  

public:
    Variant() : m_uiSelect(0) {

        BOOST_MPL_ASSERT((boost::mpl::is_sequence<Tlist>));
    };

    Variant(const Variant& var)    {
        m_spSelect = var.m_spSelect;
        m_uiSelect = var.m_uiSelect;
    }


    bool operator<(const Variant& var) const
    {
        if (m_uiSelect != var.m_uiSelect)
        {
            return (m_uiSelect < var.m_uiSelect);
        }
        else
        {
            return (m_spSelect < var.m_spSelect);
        }
    }

    bool operator==(const Variant& var) const
    {
        return (m_uiSelect == var.m_uiSelect && m_spSelect == var.m_spSelect);
    }

    bool operator!=(const Variant& var) const
    {
        return !(*this == var);
    }


protected:

    template< typename T, typename T2>
    void assign(T value) {

        typedef typename boost::mpl::find<Tlist, T2>::type index;
        typedef typename boost::mpl::end<Tlist>::type     end_index;
        BOOST_MPL_ASSERT(( boost::mpl::not_< boost::is_same<index, end_index> > ));

	std::cout << "Value in geometry" << std::endl<<value->m_storage<<std::endl;
	
        m_spSelect = boost::static_pointer_cast<void>(value);
        m_uiSelect = index::pos::value;

    };
    
    template<typename WhichType, typename CallbackType, typename VisitorType>
    typename VisitorType::result_type apply_visitor(VisitorType& rVisitor, boost::mpl::false_ /*unrolled*/) const
    {
      BOOST_ASSERT(false); //Should never assert here; only meant to stop recursion at the end of the typelist
      return typename VisitorType::result_type();
    };

  
    template<typename WhichType, typename CallbackType, typename VisitorType>
    typename VisitorType::result_type apply_visitor(VisitorType& rVisitor, boost::mpl::true_ /*unrolled*/) const
    {
      CallbackType t;
      VisitorType v;
      switch (m_uiSelect)
      {
/*#ifdef BOOST_PP_LOCAL_ITERATE
#define BOOST_PP_LOCAL_MACRO(n) \*/
      case (WhichType::value  /*+ n*/): /*\*/
        return  t.operator()<boost::mpl::int_<WhichType::value /*+ n*/> >(rVisitor,typename boost::mpl::less<boost::mpl::int_<WhichType::value /*+ n*/>, typename boost::mpl::size<Tlist>::type>::type()); /*\*/
        break;/*
#define BOOST_PP_LOCAL_LIMITS (0, 10)
#include BOOST_PP_LOCAL_ITERATE()
#endif  //BOOST_PP_LOCAL_ITERATE*/
      default:
        typedef typename boost::mpl::int_<WhichType::value + 10> next_which_t;
        return apply_visitor<next_which_t, CallbackType> ( rVisitor, typename boost::mpl::less< next_which_t, typename boost::mpl::size<Tlist>::type >::type());
      }
    };
    
    boost::shared_ptr<void> m_spSelect;
    std::size_t m_uiSelect;
};



//Unary visitation
template<typename Visitor, typename Visitable>
inline typename Visitor::result_type apply_visitor(const Visitor& visitor, Visitable& visitable)
{
    return visitable.apply_visitor<boost::mpl::int_<0> >(visitor, boost::mpl::true_());
}

template<typename Visitor, typename Visitable>
inline typename Visitor::result_type apply_visitor(Visitor& visitor, Visitable& visitable)
{
    return visitable.apply_visitor<boost::mpl::int_<0> >(visitor, boost::mpl::true_());
}

//Binary visitation
template<typename Visitor, typename Value1>
class CBinaryUnwrap2
{
public:
    typedef typename Visitor::result_type result_type;

private:
    Visitor& visitor_;
    Value1& value1_;

public:
    CBinaryUnwrap2(Visitor& visitor, Value1& value1)
            : visitor_(visitor)
            , value1_(value1)
    {
    }

public:
    template <typename Value2>
    result_type operator()(Value2& value2)
    {
        return visitor_(value1_, value2);
    }
};

template <typename Visitor, typename Visitable2>
class CBinaryUnwrap1
{
public:
    typedef typename Visitor::result_type result_type;

private:
    Visitor& visitor_;
    Visitable2& visitable2_;

public:
    CBinaryUnwrap1(Visitor& visitor, Visitable2& visitable2)
            : visitor_(visitor)
            , visitable2_(visitable2)
    {
    }

public:
    template<typename Value1>
    result_type operator()(Value1& value1)
    {
        CBinaryUnwrap2<Visitor, Value1> unwrapper(visitor_, value1);
        return apply_visitor(unwrapper, visitable2_);
    }
};

template<typename Visitor, typename Visitable1, typename Visitable2>
inline typename Visitor::result_type apply_visitor(const Visitor& visitor, Visitable1& visitable1, Visitable2& visitable2)
{
    //in 3 stappen doen; eerst visitable1 en visitable2 meegeven en daarna visit op visitable2 doen en resultaat v/d eerste meegeven, dan visit op de laatste
    dcm::CBinaryUnwrap1<const Visitor, Visitable2> unwrapper(visitor, visitable2);
    return apply_visitor(unwrapper, visitable1);
}

template<typename Visitor, typename Visitable1, typename Visitable2>
inline typename Visitor::result_type apply_visitor(Visitor& visitor, Visitable1& visitable1, Visitable2& visitable2)
{
    dcm::CBinaryUnwrap1<Visitor, Visitable2> unwrapper(visitor, visitable2);
    return apply_visitor(unwrapper, visitable1);
}

//Base class for visitor classes
template<typename R = void>
class static_visitor
{
public:
    typedef R result_type;

protected:
    // for use as base class only
    static_visitor() { }
    ~static_visitor() { }
};

template<typename Base, typename Derived>
struct is_base_of_smartptr
            : boost::is_base_of<typename Base::value_type, typename Derived::value_type>
{
};



//*****************************************
//get
//*****************************************

class failed_get : public std::exception	{

public: // std::exception implementation

    virtual const char * what() const throw()
    {
        return "get failed: "
               "the variant holds not the type requested by get<>";
    }

};

}


#endif //NS2_VARIANT_H
