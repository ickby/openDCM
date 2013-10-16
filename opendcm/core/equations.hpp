/*
    openDCM, dimensional constraint manager
    Copyright (C) 2012  Stefan Troeger <stefantroeger@gmx.net>

    This library is free software; you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more detemplate tails.

    You should have received a copy of the GNU Lesser General Public License along
    with this library; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/

#ifndef GCM_EQUATIONS_H
#define GCM_EQUATIONS_H

#include <assert.h>

#include <boost/utility/enable_if.hpp>
#include <boost/type_traits.hpp>
#include <boost/mpl/is_sequence.hpp>
#include <boost/mpl/and.hpp>
#include <boost/mpl/not.hpp>
#include <boost/fusion/include/vector.hpp>
#include <boost/fusion/include/mpl.hpp>
#include <boost/fusion/include/iterator_range.hpp>
#include <boost/fusion/include/copy.hpp>
#include <boost/fusion/include/advance.hpp>
#include <boost/fusion/include/back.hpp>
#include <boost/fusion/include/iterator_range.hpp>
#include <boost/fusion/include/nview.hpp>
#include <boost/fusion/include/for_each.hpp>
#include <boost/fusion/include/map.hpp>
#include <boost/fusion/include/as_map.hpp>

#include <boost/exception/exception.hpp>

namespace fusion = boost::fusion;
namespace mpl = boost::mpl;

#include "kernel.hpp"

namespace dcm {

struct no_option {};

template<typename Kernel>
struct Pseudo {
    typedef std::vector<typename Kernel::Vector3, Eigen::aligned_allocator<typename Kernel::Vector3> > Vec;

    template <typename DerivedA,typename DerivedB>
    void calculatePseudo(const E::MatrixBase<DerivedA>& param1, Vec& v1, const E::MatrixBase<DerivedB>& param2, Vec& v2) {};
};

template<typename Kernel>
struct Scale {
    void setScale(typename Kernel::number_type scale) {};
};

template<typename Kernel>
struct PseudoScale {
    typedef std::vector<typename Kernel::Vector3, Eigen::aligned_allocator<typename Kernel::Vector3> > Vec;

    template <typename DerivedA,typename DerivedB>
    void calculatePseudo(const E::MatrixBase<DerivedA>& param1, Vec& v1, const E::MatrixBase<DerivedB>& param2, Vec& v2) {};
    void setScale(typename Kernel::number_type scale) {};
};

//type to allow a metaprogramming check for a Equation
struct EQ {};

template<typename Seq, typename T>
struct pushed_seq;

template<typename seq>
struct constraint_sequence : public seq {

    //an equation gets added to this equation
    template<typename T>
    typename boost::enable_if< boost::is_base_of< dcm::EQ, T>, typename pushed_seq<seq, T>::type >::type operator &(T val) {

        typedef typename pushed_seq<seq, T>::type Sequence;
        typedef typename fusion::result_of::begin<Sequence>::type Begin;
        typedef typename fusion::result_of::find<Sequence, typename fusion::result_of::back<typename pushed_seq<seq, T>::S1>::type >::type EndOld;


        //create the new sequence
        Sequence vec;

        //copy the old values into the new sequence
        Begin b(vec);
        EndOld eo(vec);

        fusion::iterator_range<Begin, EndOld> range(b, eo);
        fusion::copy(*this, range);

        //insert this object at the end of the sequence
        *fusion::find<T>(vec) = val;

        //and return our new extendet sequence
        return vec;
    };

    template<typename T>
    void pretty(T type) {
        std::cout<<"pretty: "<<__PRETTY_FUNCTION__<<std::endl;
    };

    //an sequence gets added to this equation (happens only if sequenced equations like coincident are used)
    template<typename T>
    typename boost::enable_if< mpl::is_sequence<T>, typename pushed_seq<T, seq>::type >::type operator &(T val) {

        typedef typename pushed_seq<T, seq>::type Sequence;
        typedef typename fusion::result_of::begin<Sequence>::type Begin;
        typedef typename fusion::result_of::find<Sequence, typename fusion::result_of::back<typename pushed_seq<T, seq>::S1>::type >::type EndF;


        //create the new sequence
        Sequence vec;

        Begin b(vec);
        EndF ef(vec);

        fusion::iterator_range<Begin, EndF> range(b, ef);
        fusion::copy(val, range);


        //to copy the types of the second sequence is not as easy as before. If types were already present in
        //the original sequence they are not added again. therefore we need to find all types of the second sequence
        //in the new one and assign the objects to this positions.

        //get a index vector for all second-sequence-elements
        typedef typename mpl::transform<typename pushed_seq<T, seq>::S2,
                fusion::result_of::distance<typename fusion::result_of::begin<Sequence>::type,
                fusion::result_of::find<Sequence, mpl::_1> > >::type position_vector;

        //and copy the types in
        fusion::nview<Sequence, position_vector> view(vec);
        fusion::copy(*this, view);

        //and return our new extendet sequence
        return vec;
    };
};

template<typename Seq, typename T>
struct pushed_seq {
    typedef typename mpl::if_<mpl::is_sequence<Seq>, Seq, fusion::vector1<Seq> >::type S1;
    typedef typename mpl::if_<mpl::is_sequence<T>, T, fusion::vector1<T> >::type S2;

    typedef typename mpl::fold< S2, S1, mpl::if_< boost::is_same<
    mpl::find<mpl::_1, mpl::_2>, mpl::end<mpl::_1> >, mpl::push_back<mpl::_1,mpl::_2>, mpl::_1> >::type unique_vector;

    typedef typename fusion::result_of::as_vector< unique_vector >::type vec;
    typedef constraint_sequence<vec> type;
};

template<typename Derived, typename Option, bool rotation_only = false>
struct Equation : public EQ {

    typedef typename mpl::if_<mpl::is_sequence<Option>, Option, mpl::vector<Option> >::type option_sequence;
    typedef typename mpl::fold<option_sequence, fusion::map<>, fusion::result_of::push_back<mpl::_1, fusion::pair<mpl::_2, std::pair<bool, mpl::_2> > > > ::type option_set_map;
    typedef typename fusion::result_of::as_map<option_set_map>::type options;

    options values;
    bool pure_rotation;

    struct option_copy {

        options& values;
        option_copy(options& op) : values(op) {};

        template<typename T>
        void operator()(const T& val) {
            if(val.first)
                fusion::at_key<T>(values) = val;
        };
    };

    Equation() : pure_rotation(rotation_only) {};

    template<typename T>
    typename boost::enable_if<fusion::result_of::has_key<options, T>, Derived&>::type operator()(const T& val) {
        fusion::at_key<T>(values).second = val;
        fusion::at_key<T>(values).first  = true;
        return *(static_cast<Derived*>(this));
    };
    //assign option
    template<typename T>
    typename boost::enable_if<fusion::result_of::has_key<options, T>, Derived&>::type operator=(const T& val) {
        return operator()(val);
    };
    //assign complete equation
    template<typename T>
    typename boost::enable_if<boost::is_base_of<EQ, T>, Derived& >::type
    operator=(const T& eq) {

        //we only copy the values which were set and are therefore valid
	option_copy oc(values);
	fusion::for_each(eq.values, oc);
    };

    //an equation gets added to this equation
    template<typename T>
    typename boost::enable_if< boost::is_base_of< dcm::EQ, T>, typename pushed_seq<T, Derived>::type >::type operator &(T val) {

        typename pushed_seq<T, Derived>::type vec;
        *fusion::find<T>(vec) = val;
        *fusion::find<Derived>(vec) = *(static_cast<Derived*>(this));
        return vec;
    };

    //an sequence gets added to this equation (happens only if sequenced equations like coincident are used)
    template<typename T>
    typename boost::enable_if< mpl::is_sequence<T>, typename pushed_seq<T, Derived>::type >::type operator &(T val) {

        typedef typename pushed_seq<T, Derived>::type Sequence;
        typedef typename fusion::result_of::begin<Sequence>::type Begin;
        typedef typename fusion::result_of::find<Sequence, typename fusion::result_of::back<typename pushed_seq<T, Derived>::S1>::type >::type EndOld;

        //create the new sequence
        Sequence vec;

        //copy the old values into the new sequence
        Begin b(vec);
        EndOld eo(vec);

        fusion::iterator_range<Begin, EndOld> range(b, eo);
        fusion::copy(val, range);

        //insert this object at the end of the sequence
        *fusion::find<Derived>(vec) = *static_cast<Derived*>(this);

        //and return our new extendet sequence
        return vec;
    };
};

struct Distance : public Equation<Distance, double> {

    using Equation::operator=;
    using Equation::options;
    Distance() : Equation() {};

    template< typename Kernel, typename Tag1, typename Tag2 >
    struct type {

        type() {
            throw constraint_error() << boost::errinfo_errno(100) << error_message("unsupported geometry in distance constraint")
                                     << error_type_first_geometry(typeid(Tag1).name()) << error_type_second_geometry(typeid(Tag2).name());
        };

        typedef typename Kernel::number_type Scalar;
        typedef typename Kernel::VectorMap   Vector;
        typedef std::vector<typename Kernel::Vector3, Eigen::aligned_allocator<typename Kernel::Vector3> > Vec;

        options values;
        //template definition
        template <typename DerivedA,typename DerivedB>
        void calculatePseudo(const E::MatrixBase<DerivedA>& param1, Vec& v1, const E::MatrixBase<DerivedB>& param2, Vec& v2) {
            assert(false);
        };
        void setScale(Scalar scale) {
            assert(false);
        };
        template <typename DerivedA,typename DerivedB>
        Scalar calculate(const E::MatrixBase<DerivedA>& param1,  const E::MatrixBase<DerivedB>& param2) {
            assert(false);
            return 0;
        };
        template <typename DerivedA,typename DerivedB, typename DerivedC>
        Scalar calculateGradientFirst(const E::MatrixBase<DerivedA>& param1,
                                      const E::MatrixBase<DerivedB>& param2,
                                      const E::MatrixBase<DerivedC>& dparam1) {
            assert(false);
            return 0;
        };
        template <typename DerivedA,typename DerivedB, typename DerivedC>
        Scalar calculateGradientSecond(const E::MatrixBase<DerivedA>& param1,
                                       const E::MatrixBase<DerivedB>& param2,
                                       const E::MatrixBase<DerivedC>& dparam2) {
            assert(false);
            return 0;
        };
        template <typename DerivedA,typename DerivedB, typename DerivedC>
        void calculateGradientFirstComplete(const E::MatrixBase<DerivedA>& param1,
                                            const E::MatrixBase<DerivedB>& param2,
                                            E::MatrixBase<DerivedC>& gradient) {
            assert(false);
        };
        template <typename DerivedA,typename DerivedB, typename DerivedC>
        void calculateGradientSecondComplete(const E::MatrixBase<DerivedA>& param1,
                                             const E::MatrixBase<DerivedB>& param2,
                                             E::MatrixBase<DerivedC>& gradient) {
            assert(false);
        };
    };
};

//the possible directions
enum Direction { parallel, equal, opposite, perpendicular };

struct Orientation : public Equation<Orientation, Direction, true> {

    using Equation::operator=;
    using Equation::options;
    Orientation() : Equation() {};

    template< typename Kernel, typename Tag1, typename Tag2 >
    struct type : public PseudoScale<Kernel> {

        type() {
            throw constraint_error() << boost::errinfo_errno(101) << error_message("unsupported geometry in orientation constraint")
                                     << error_type_first_geometry(typeid(Tag1).name()) << error_type_second_geometry(typeid(Tag2).name());
        };

        typedef typename Kernel::number_type Scalar;
        typedef typename Kernel::VectorMap   Vector;

        options values;

        //template definition
        template <typename DerivedA,typename DerivedB>
        Scalar calculate(const E::MatrixBase<DerivedA>& param1,  const E::MatrixBase<DerivedB>& param2) {
            assert(false);
            return 0;
        };
        template <typename DerivedA,typename DerivedB, typename DerivedC>
        Scalar calculateGradientFirst(const E::MatrixBase<DerivedA>& param1,
                                      const E::MatrixBase<DerivedB>& param2,
                                      const E::MatrixBase<DerivedC>& dparam1) {
            assert(false);
            return 0;
        };
        template <typename DerivedA,typename DerivedB, typename DerivedC>
        Scalar calculateGradientSecond(const E::MatrixBase<DerivedA>& param1,
                                       const E::MatrixBase<DerivedB>& param2,
                                       const E::MatrixBase<DerivedC>& dparam2) {
            assert(false);
            return 0;
        };
        template <typename DerivedA,typename DerivedB, typename DerivedC>
        void calculateGradientFirstComplete(const E::MatrixBase<DerivedA>& param1,
                                            const E::MatrixBase<DerivedB>& param2,
                                            E::MatrixBase<DerivedC>& gradient) {
            assert(false);
        };
        template <typename DerivedA,typename DerivedB, typename DerivedC>
        void calculateGradientSecondComplete(const E::MatrixBase<DerivedA>& param1,
                                             const E::MatrixBase<DerivedB>& param2,
                                             E::MatrixBase<DerivedC>& gradient) {
            assert(false);
        };
    };
};

struct Angle : public Equation<Angle, double, true> {

    using Equation::operator=;
    Angle() : Equation() {};

    template< typename Kernel, typename Tag1, typename Tag2 >
    struct type : public PseudoScale<Kernel> {

        type() {
            throw constraint_error() << boost::errinfo_errno(102) << error_message("unsupported geometry in angle constraint")
                                     << error_type_first_geometry(typeid(Tag1).name()) << error_type_second_geometry(typeid(Tag2).name());
        };

        typedef typename Kernel::number_type Scalar;
        typedef typename Kernel::VectorMap   Vector;

        options values;

        //template definition
        template <typename DerivedA,typename DerivedB>
        Scalar calculate(const E::MatrixBase<DerivedA>& param1,  const E::MatrixBase<DerivedB>& param2) {
            assert(false);
            return 0;
        };
        template <typename DerivedA,typename DerivedB, typename DerivedC>
        Scalar calculateGradientFirst(const E::MatrixBase<DerivedA>& param1,
                                      const E::MatrixBase<DerivedB>& param2,
                                      const E::MatrixBase<DerivedC>& dparam1) {
            assert(false);
            return 0;
        };
        template <typename DerivedA,typename DerivedB, typename DerivedC>
        Scalar calculateGradientSecond(const E::MatrixBase<DerivedA>& param1,
                                       const E::MatrixBase<DerivedB>& param2,
                                       const E::MatrixBase<DerivedC>& dparam2) {
            assert(false);
            return 0;
        };
        template <typename DerivedA,typename DerivedB, typename DerivedC>
        void calculateGradientFirstComplete(const E::MatrixBase<DerivedA>& param1,
                                            const E::MatrixBase<DerivedB>& param2,
                                            E::MatrixBase<DerivedC>& gradient) {
            assert(false);
        };
        template <typename DerivedA,typename DerivedB, typename DerivedC>
        void calculateGradientSecondComplete(const E::MatrixBase<DerivedA>& param1,
                                             const E::MatrixBase<DerivedB>& param2,
                                             E::MatrixBase<DerivedC>& gradient) {
            assert(false);
        };
    };
};

//static is needed to restrain the scope of the objects to the current compilation unit. Without it
//every compiled file including this header would define these as global and the linker would find
//multiple definitions of the same objects
static Distance distance;
static Orientation orientation;
static Angle    angle;

};

#endif //GCM_EQUATIONS_H

