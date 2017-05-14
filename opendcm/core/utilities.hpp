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
struct SequenceApplyer {
    
    Functor<Sequence> functor;
    
    template<typename T>
    SequenceApplyer(T& param) 
        : functor(Functor<Sequence>(param)) {};
        
    SequenceApplyer(SequenceApplyer& r) 
        : functor(r.functor) {};
        
    template<typename T>
    void operator()(const T& t) {        
        functor.template operator()<T>();
    };
};

template<typename Sequence, typename Functor>
struct RecursiveSequenceApplyer {
    
    Functor functor;
    
    template<typename... Args>
    RecursiveSequenceApplyer(Args&... param) 
        : functor(Functor(param...)) {};
    
    RecursiveSequenceApplyer(RecursiveSequenceApplyer& r) 
        : functor(r.functor) {};
        
    template<typename T>
    void operator()(const T& t) {        
            boost::mpl::for_each<Sequence>(InnerLoop<T>(functor));
    };
    
    template<typename T1>
    struct InnerLoop {

        Functor& functor;
        
        InnerLoop(Functor& f) 
            : functor(f) {};
            
        template<typename T2>
        void operator()(const T2& t) {            
            functor.template operator()<T1, T2>();
        };
    };
};

class IDGeneratorBase {
protected:      
    static int generateId(){
        static int idx = -1;
        return ++idx;
    }
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







