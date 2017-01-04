/*
    openDCM, dimensional constraint manager
    Copyright (C) 2014  Stefan Troeger <stefantroeger@gmx.net>

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

#ifndef DCM_UTILITIES_H
#define DCM_UTILITIES_H

#include <boost/variant.hpp>
#include <boost/mpl/for_each.hpp>
#include <boost/mpl/if.hpp>

namespace dcm {
namespace utilities {
    
template<typename ...types>
struct Variant {
    
    //Create the variant. Note that we want to make it possible to detect if the variant
    //holds nothing, for this purpose we add a boost::blank in front of the typelist. The 
    //first type is the one set by default from the variant.
    typedef boost::variant< boost::blank, types... > VariantType;
    
    template<typename Visitor>
    typename Visitor::result_type apply(Visitor& vis) {
        return boost::apply_visitor(vis, m_variant);
    };
    template<typename Visitor>
    typename Visitor::result_type apply(const Visitor& vis) const{
        return boost::apply_visitor(vis, m_variant);
    };
    
    template<typename T>
    typename boost::add_reference<T>::type get() {
        return boost::get<T>(m_variant);
    };
        
    bool holdsType() {
        return m_variant.which()!=0;
    };
    
protected:
    VariantType m_variant;
    
    void clearVariant() {
        m_variant = boost::blank();
    };
};

template<typename Sequence, template<class> class Functor>
struct RecursiveSequenceApplyer {
    
    Functor<Sequence> functor;
    
    template<typename T>
    RecursiveSequenceApplyer(T& param) 
        : functor(Functor<Sequence>(param)) {};
    
    RecursiveSequenceApplyer(RecursiveSequenceApplyer& r) 
        : functor(r.functor) {};
        
    template<typename T>
    void operator()(const T& t) {        
            typedef mpl::range_c<int, T::value, mpl::size<Sequence>::value> StorageRange;
            boost::mpl::for_each<StorageRange>(InnerLoop<T>(functor));
    };
    
    template<typename Number>
    struct InnerLoop {

        Functor<Sequence>& functor;
        
        InnerLoop(Functor<Sequence>& f) 
            : functor(f) {};
            
        template<typename T>
        void operator()(const T& t) {            
            functor.template operator()<Number, T>();
        };
    };
};

//generates a int index value of a type in a sequence
template<typename Sequence, typename G>
struct index {
    typedef typename mpl::find<Sequence, G>::type iterator;
    typedef boost::is_same<iterator, typename mpl::end<Sequence>::type> valid;
    BOOST_MPL_ASSERT_MSG(mpl::not_<valid>::value, CONSTRAINT_TYPE_NOT_REGISTERT, (G));
    
    typedef typename mpl::distance<typename mpl::begin<Sequence>::type,
                        iterator>::type type; 
    const static long value = type::value;
};

//exposes the given type except if it is floating point, than int is exposed
template<typename T>
struct non_floating : public boost::mpl::if_<std::is_floating_point<T>, int, T>{};
    
}//utilities

namespace details {
//allow direct access to the stored geometry in any TypeVariant, copyed from boost variant get
    
template <typename T>
struct get_visitor {
private:

    typedef typename boost::add_pointer<T>::type pointer;
    typedef typename boost::add_reference<T>::type reference;

public:
    typedef pointer result_type;

public:
    pointer operator()(reference operand) const   {
        return boost::addressof(operand);
    }

    template <typename U>
    pointer operator()(const U&) const  {
        return static_cast<pointer>(0);
    }
};
} //details

template<typename T, typename Types>
typename boost::add_reference<T>::type 
get(std::shared_ptr<utilities::Variant<Types>> variant) {

    typedef typename boost::add_pointer<T>::type T_ptr;
    details::get_visitor<T> v;
    T_ptr result = variant->apply(v);

    if (!result)
        throw boost::bad_get();
    return *result;
};

template<typename T = void>
using visitor = boost::static_visitor<T>;

}//dcm

#endif //DCM_UTILITIES_H







