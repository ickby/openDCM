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
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License along
    with this library; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/

#ifndef DCM_TEST_MODULE3D_HPP
#define DCM_TEST_MODULE3D_HPP

#include "opendcm/core.hpp"
#include "opendcm/module3d.hpp"

#include "test/Octave/debugsolver.hpp"

#include <time.h>
#include <iostream>
#include <iomanip>

#include <boost/test/unit_test.hpp>

struct point : std::vector<double> {};
typedef Eigen::Matrix<double, 6,1> line_t;

namespace dcm {

template<>
struct geometry_traits<point> {
    typedef tag::direction3D tag;
    typedef modell::XYZ modell;
    typedef orderd_bracket_accessor accessor;
};

template<>
struct geometry_traits<Eigen::Vector3d> {
    typedef tag::point3D tag;
    typedef modell::XYZ modell;
    typedef orderd_roundbracket_accessor accessor;
};


template<>
struct geometry_traits<line_t> {
    typedef tag::line3D  tag;
    typedef modell::XYZ2 modell;
    typedef orderd_roundbracket_accessor accessor;
};

}

//two vectors perpendicular, maybe the easiest constraints of them all
struct test_constraint : public dcm::Equation<test_constraint, int, 99> {

    using Equation::options;
    using Equation::values;

    void setDefault() {
        fusion::at_key<int>(values) = std::make_pair(false, 0);
    };

    template< typename Kernel, typename Tag1, typename Tag2 >
    struct type : public dcm::PseudoScale<Kernel> {
        typedef typename Kernel::number_type Scalar;
        typedef typename Kernel::VectorMap   Vector;
        typename test_constraint::options values;

        template <typename DerivedA,typename DerivedB>
        Scalar calculate(const E::MatrixBase<DerivedA>& param1,  const E::MatrixBase<DerivedB>& param2)  {
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
                                       const E::MatrixBase<DerivedC>& dparam2)  {
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


    template< typename Kernel >
    struct type<Kernel, dcm::tag::direction3D, dcm::tag::direction3D> : public dcm::PseudoScale<Kernel> {

        typedef typename Kernel::number_type Scalar;
        typedef typename Kernel::VectorMap   Vector;
        typename test_constraint::options values;

        template <typename DerivedA,typename DerivedB>
        Scalar calculate(const E::MatrixBase<DerivedA>& param1,  const E::MatrixBase<DerivedB>& param2) {
            return param1.dot(param2);
        };

        template <typename DerivedA,typename DerivedB, typename DerivedC>
        Scalar calculateGradientFirst(const E::MatrixBase<DerivedA>& param1,
                                      const E::MatrixBase<DerivedB>& param2,
                                      const E::MatrixBase<DerivedC>& dparam1) {

            return dparam1.dot(param2);
        };

        template <typename DerivedA,typename DerivedB, typename DerivedC>
        Scalar calculateGradientSecond(const E::MatrixBase<DerivedA>& param1,
                                       const E::MatrixBase<DerivedB>& param2,
                                       const E::MatrixBase<DerivedC>& dparam2) {

            return param1.dot(dparam2);
        };

        template <typename DerivedA,typename DerivedB, typename DerivedC>
        void calculateGradientFirstComplete(const E::MatrixBase<DerivedA>& param1,
                                            const E::MatrixBase<DerivedB>& param2,
                                            E::MatrixBase<DerivedC>& gradient) {

            gradient(0) = param2(0);
            gradient(1) = param2(1);
            gradient(2) = param2(2);
        };

        template <typename DerivedA,typename DerivedB, typename DerivedC>
        void calculateGradientSecondComplete(const E::MatrixBase<DerivedA>& param1,
                                             const E::MatrixBase<DerivedB>& param2,
                                             E::MatrixBase<DerivedC>& gradient) {

            gradient(0) = param1(0);
            gradient(1) = param1(1);
            gradient(2) = param1(2);
        };
    };

    template<typename Kernel>
    struct type<Kernel, dcm::tag::point3D, dcm::tag::point3D> : public type<Kernel, dcm::tag::direction3D, dcm::tag::direction3D> {};

};

//multi-equation constraint test
typedef fusion::vector2<test_constraint, dcm::Distance> vector;
struct comp_constraint : public dcm::constraint_sequence<vector> {
    //allow to set the distance
    comp_constraint& operator()(double val) {
        fusion::at_c<1>(*this) = val;
        return *this;
    };
};

typedef dcm::Kernel<double> kernel;
typedef dcm::Module3D< mpl::vector3<point, Eigen::Vector3d, line_t > > Module;
typedef dcm::Module3D< mpl::vector3<point, Eigen::Vector3d, line_t >, std::string > ModuleID;
typedef dcm::System<kernel, Module> SystemNOID;
typedef dcm::System<kernel, ModuleID> SystemID;
typedef Module::type<SystemNOID>::Geometry3D geom;
typedef ModuleID::type<SystemID>::Geometry3D geomid;
typedef boost::shared_ptr<geom> geom_ptr;
typedef boost::shared_ptr<geomid> geomid_ptr;

typedef Module::type<SystemNOID>::Constraint3D cons;
typedef ModuleID::type<SystemID>::Constraint3D consid;
typedef boost::shared_ptr<cons> cons_ptr;
typedef boost::shared_ptr<consid> consid_ptr;

typedef SystemNOID::Cluster::vertex_iterator viter;
typedef Module::type<SystemNOID>::vertex_prop vertex_prop;

#endif