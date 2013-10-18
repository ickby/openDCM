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

#include "opendcm/core.hpp"
#include "opendcm/module3d.hpp"
#include "opendcm/modulepart.hpp"
#include "opendcm/moduleshape3d.hpp"

#include <boost/mpl/vector.hpp>
#include <boost/test/unit_test.hpp>

namespace mpl = boost::mpl;

typedef Eigen::Matrix<double, 3,1> point_t;
struct line_t : public Eigen::Matrix<double, 6,1> {};
struct plane_t : public Eigen::Matrix<double, 6,1> {};
struct segment_t : public Eigen::Matrix<double, 6,1> {};
struct cylinder_t : public Eigen::Matrix<double, 7,1> {};

struct place {
    Eigen::Quaterniond quat;
    Eigen::Vector3d trans;

    place():quat(1,2,3,4) {
        quat.normalize();
        trans<<1,2,3;
    };
};
struct place_accessor {

    template<typename Scalar, int ID, typename T>
    Scalar get(T& t) {
        switch(ID) {
        case 0:
            return t.quat.w();
        case 1:
            return t.quat.x();
        case 2:
            return t.quat.y();
        case 3:
            return t.quat.z();
        case 4:
            return t.trans(0);
        case 5:
            return t.trans(1);
        case 6:
            return t.trans(2);
        default:
            return 0;
        };
    };
    template<typename Scalar, int ID, typename T>
    void set(Scalar value, T& t) {
        switch(ID) {
        case 0:
            t.quat.w() = value;
            break;
        case 1:
            t.quat.x() = value;
            break;
        case 2:
            t.quat.y() = value;
            break;
        case 3:
            t.quat.z() = value;
            break;
        case 4:
            t.trans(0) = value;
            break;
        case 5:
            t.trans(1) = value;
            break;
        case 6:
            t.trans(2) = value;
            break;
        };
    };
    template<typename T>
    void finalize(T& t) {};
};


namespace dcm {

template<>
struct geometry_traits<point_t> {
    typedef tag::point3D tag;
    typedef modell::XYZ modell;
    typedef orderd_bracket_accessor accessor;
};

template<>
struct geometry_traits<plane_t> {
    typedef tag::plane3D tag;
    typedef modell::XYZ2 modell;
    typedef orderd_roundbracket_accessor accessor;
};


template<>
struct geometry_traits<line_t> {
    typedef tag::line3D  tag;
    typedef modell::XYZ2 modell;
    typedef orderd_roundbracket_accessor accessor;
};

template<>
struct geometry_traits<segment_t> {
    typedef tag::segment3D  tag;
    typedef modell::XYZ2 modell;
    typedef orderd_roundbracket_accessor accessor;
};

template<>
struct geometry_traits<cylinder_t> {
    typedef tag::cylinder3D  tag;
    typedef modell::XYZ2P modell;
    typedef orderd_roundbracket_accessor accessor;
};

template<>
struct geometry_traits< place > {
    typedef tag::part  tag;
    typedef modell::quaternion_wxyz_vec3 modell;
    typedef place_accessor accessor;
};

}

//solver for pure differential checking
template<typename Kernel>
struct CheckSolver {

    typedef typename Kernel::number_type Scalar;
    CheckSolver(Kernel* k) {};

    template<typename Functor>
    int solve(typename Kernel::MappedEquationSystem& sys, Functor& f) {

        int jcount = 100;
        double diffmax = 0;
        int fails = 0;
        for(int i=0; i<sys.Parameter.rows(); i++) {

            //varying parameter i
            sys.Parameter(i) -= 0.1*double(jcount/2);
            for(int j=0; j<jcount; j++) {

                //get first new value
                sys.Parameter(i) += 0.1 - 1e-5;
                sys.recalculate();
                double dleft = sys.Residual(0);

                //get second new value
                sys.Parameter(i) += 2e-5;
                sys.recalculate();
                double dright = sys.Residual(0);

                //get the residual diff between the last two calculations
                double rdiff = (dright - dleft)/2e-5;

                diffmax = std::max(sys.Jacobi(0,i), diffmax);

                //compare with the algeraic caclculated diff
                if(std::abs(sys.Jacobi(0,i) - rdiff) > 1e-3) {
                    std::cout<<"iter: "<<j<<", parameter: "<<i<<", calc: "<<sys.Jacobi(0,i)<<" vs real: "<<rdiff<<std::endl;
                    fails++;
                    if(fails>3)
                        throw(std::exception()); //and throw if the difference is too much
                }
                //set the parameter to the basic value to avoid drift
                sys.Parameter(i) -= 1e-5;
            };
            //set the initial value
            sys.Parameter(i) -= 0.1*double(jcount/2);
        };
        return 1;
    }
};



using namespace dcm;

BOOST_AUTO_TEST_SUITE(constraint3d_test_suit);

typedef dcm::Kernel<double, CheckSolver> Kernel;
typedef Module3D< mpl::vector5<point_t, line_t, plane_t, cylinder_t, segment_t > > Module;
typedef ModuleShape3D< mpl::vector0<> > ModuleShape; //need that to allow shapes
typedef dcm::ModulePart< mpl::vector1< place > > ModulePart;
typedef System<Kernel, Module, ModuleShape, ModulePart> System;

typedef Module::type<System>::Geometry3D geom;
typedef boost::shared_ptr<geom> geom_ptr;

typedef Module::type<System>::Constraint3D cons;
typedef boost::shared_ptr<cons> cons_ptr;

typedef ModulePart::type<System>::Part part;
typedef boost::shared_ptr<part> part_ptr;

typedef Module::type<System>::vertex_prop vertex_prop;
typedef System::Cluster ClusterGraph;

struct notype {};
template<typename t1, typename t2, typename ct>
struct createConstraint {
    cons_ptr apply(System& system, geom_ptr g1, geom_ptr g2, t1 val1, t2 val2, ct type) {
        return system.createConstraint3D(g2, g1, (type=val1) & (type=val2));
    };
};

template<typename t1, typename ct>
struct createConstraint<t1, notype, ct> {
    cons_ptr apply(System& system, geom_ptr g1, geom_ptr g2, t1 val1, notype val2, ct type) {
        return system.createConstraint3D(g1, g2, type=val1);
    };
};

//automated checking of differentials
template<typename T1, typename T2, typename CT>
struct constraint_checker {

    System& system;

    constraint_checker(System& sys) : system(sys) { };

    template<typename TT1, typename TT2>
    bool check_normal(TT1 val1, TT2 val2 = TT2()) {

        //create the geometries
        geom_ptr g1, g2;
        T1 t1;
        t1.setRandom();
        T2 t2;
        t2.setRandom();
        g1 = system.createGeometry3D(t1);
        g2 = system.createGeometry3D(t2);

        //create the constraint and set the wanted value
        CT type;
        createConstraint<TT1, TT2, CT>().apply(system,g1,g2,val1,val2,type);

        try {
            system.solve();
            system.removeGeometry3D(g1);
            system.removeGeometry3D(g2);
            return true;
        }
        catch(...) {
            system.removeGeometry3D(g1);
            system.removeGeometry3D(g2);
            return false;
        }
    }

    template<typename TT1, typename TT2>
    bool check_cluster(TT1 val1, TT2 val2 = TT2()) {

        place pl1;
        pl1.trans *= 100;
        part_ptr p1 = system.createPart(pl1);
        part_ptr p2 = system.createPart(place());

        //create the geometries
        geom_ptr g1, g2;
        T1 t1;
        t1.setRandom();
        T2 t2;
        t2.setRandom();
        g1 = p1->addGeometry3D(t1);
        g2 = p2->addGeometry3D(t2);

        //create the constraint and set the wanted value
        CT type;
        createConstraint<TT1, TT2, CT>().apply(system,g1,g2,val1,val2,type);

        try {
            system.solve();
            system.removePart(p1);
            system.removePart(p2);
            return true;
        }
        catch(...) {
            system.removePart(p1);
            system.removePart(p2);
            return false;
        }
    }

};

template<typename T1, typename T2>
struct constraint_checker_orientation  {

    System& sys;
    constraint_checker_orientation(System& v) : sys(v) { };

    bool check() {

        bool b1,b2,b3,b4,b5,b6,b7,b8;

        constraint_checker<T1, T2, dcm::Orientation> checker(sys);
        BOOST_CHECK(b1=checker.check_normal(dcm::equal, notype()));
        BOOST_CHECK(b2=checker.check_cluster(dcm::equal, notype()));

        constraint_checker<T1, T2, dcm::Orientation> checker1(sys);
        BOOST_CHECK(b3=checker1.check_normal(dcm::opposite, notype()));
        BOOST_CHECK(b4=checker1.check_cluster(dcm::opposite, notype()));

        constraint_checker<T1, T2, dcm::Orientation> checker2(sys);
        BOOST_CHECK(b5=checker2.check_normal(dcm::parallel, notype()));
        BOOST_CHECK(b6=checker2.check_cluster(dcm::parallel, notype()));

        constraint_checker<T1, T2, dcm::Orientation> checker3(sys);
        BOOST_CHECK(b7=checker3.check_normal(dcm::perpendicular, notype()));
        BOOST_CHECK(b8=checker3.check_cluster(dcm::perpendicular, notype()));

        return b1 && b2 && b3 && b4 && b5 && b6 && b7 && b8;
    };
};


BOOST_AUTO_TEST_CASE(constraint3d_distance) {

    System sys;
    constraint_checker<point_t, point_t, dcm::Distance> checker(sys);
    BOOST_CHECK(checker.check_normal(2., notype()));
    BOOST_CHECK(checker.check_cluster(2., notype()));

    constraint_checker<point_t, line_t, dcm::Distance> checker2(sys);
    BOOST_CHECK(checker2.check_normal(2., notype()));
    BOOST_CHECK(checker2.check_cluster(2., notype()));

    constraint_checker<point_t, plane_t, dcm::Distance> checker3(sys);
    BOOST_CHECK(checker3.check_normal(2., notype()));
    BOOST_CHECK(checker3.check_cluster(2., notype()));

    constraint_checker<point_t, cylinder_t, dcm::Distance> checker4(sys);
    BOOST_CHECK(checker4.check_normal(2., notype()));
    BOOST_CHECK(checker4.check_cluster(2., notype()));

    constraint_checker<line_t, line_t, dcm::Distance> checker5(sys);
    BOOST_CHECK(checker5.check_normal(2., notype()));
    BOOST_CHECK(checker5.check_cluster(2., notype()));

    constraint_checker<line_t, plane_t, dcm::Distance> checker6(sys);
    BOOST_CHECK(checker6.check_normal(2., notype()));
    BOOST_CHECK(checker6.check_cluster(2., notype()));

    constraint_checker<line_t, cylinder_t, dcm::Distance> checker7(sys);
    BOOST_CHECK(checker7.check_normal(2., notype()));
    BOOST_CHECK(checker7.check_cluster(2., notype()));

    constraint_checker<plane_t, plane_t, dcm::Distance> checker8(sys);
    BOOST_CHECK(checker8.check_normal(2., notype()));
    BOOST_CHECK(checker8.check_cluster(2., notype()));

    constraint_checker<plane_t, cylinder_t, dcm::Distance> checker9(sys);
    BOOST_CHECK(checker9.check_normal(2., notype()));
    BOOST_CHECK(checker9.check_cluster(2., notype()));

    constraint_checker<cylinder_t, cylinder_t, dcm::Distance> checker10(sys);
    BOOST_CHECK(checker10.check_normal(2., notype()));
    BOOST_CHECK(checker10.check_cluster(2., notype()));
    
    //check different solution spaces where possible
    constraint_checker<point_t, plane_t, dcm::Distance> checker11(sys);
    BOOST_CHECK(checker11.check_normal(2., dcm::positiv_directional));
    BOOST_CHECK(checker11.check_cluster(2., dcm::positiv_directional));
    
    constraint_checker<point_t, plane_t, dcm::Distance> checker12(sys);
    BOOST_CHECK(checker12.check_normal(2., dcm::negative_directional));
    BOOST_CHECK(checker12.check_cluster(2., dcm::negative_directional));
    
    constraint_checker<point_t, cylinder_t, dcm::Distance> checker13(sys);
    BOOST_CHECK(checker13.check_normal(2., dcm::positiv_directional));
    BOOST_CHECK(checker13.check_cluster(2., dcm::positiv_directional));
    
    constraint_checker<point_t, cylinder_t, dcm::Distance> checker14(sys);
    BOOST_CHECK(checker14.check_normal(2., dcm::negative_directional));
    BOOST_CHECK(checker14.check_cluster(2., dcm::negative_directional));
}
/*
BOOST_AUTO_TEST_CASE(constraint3d_orientation) {

    System sys;
    constraint_checker_orientation<line_t, line_t> checker(sys);
    BOOST_CHECK(checker.check());

    constraint_checker_orientation<line_t, plane_t> checker1(sys);
    BOOST_CHECK(checker1.check());

    constraint_checker_orientation<line_t, cylinder_t> checker2(sys);
    BOOST_CHECK(checker2.check());

    constraint_checker_orientation<plane_t, plane_t> checker3(sys);
    BOOST_CHECK(checker3.check());

    constraint_checker_orientation<plane_t, cylinder_t> checker4(sys);
    BOOST_CHECK(checker4.check());

    constraint_checker_orientation<cylinder_t, cylinder_t> checker5(sys);
    BOOST_CHECK(checker5.check());
}

BOOST_AUTO_TEST_CASE(constraint3d_angle) {

    System sys;
    constraint_checker<line_t, line_t, dcm::Angle> checker(sys);
    BOOST_CHECK(checker.check_normal(2., notype()));
    BOOST_CHECK(checker.check_cluster(2., notype()));

    constraint_checker<line_t, plane_t, dcm::Angle> checker1(sys);
    BOOST_CHECK(checker1.check_normal(2., notype()));
    BOOST_CHECK(checker1.check_cluster(2., notype()));

    constraint_checker<line_t, cylinder_t, dcm::Angle> checker2(sys);
    BOOST_CHECK(checker2.check_normal(2., notype()));
    BOOST_CHECK(checker2.check_cluster(2., notype()));

    constraint_checker<plane_t, plane_t, dcm::Angle> checker3(sys);
    BOOST_CHECK(checker3.check_normal(2., notype()));
    BOOST_CHECK(checker3.check_cluster(2., notype()));

    constraint_checker<plane_t, cylinder_t, dcm::Angle> checker4(sys);
    BOOST_CHECK(checker4.check_normal(2., notype()));
    BOOST_CHECK(checker4.check_cluster(2., notype()));

    constraint_checker<cylinder_t, cylinder_t, dcm::Angle> checker5(sys);
    BOOST_CHECK(checker5.check_normal(2., notype()));
    BOOST_CHECK(checker5.check_cluster(2., notype()));
}

BOOST_AUTO_TEST_CASE(constraint3d_shape_distance) {

    System sys;
    constraint_checker<point_t, segment_t, dcm::Distance> checker(sys);
    BOOST_CHECK(checker.check_normal(2., notype()));
    BOOST_CHECK(checker.check_cluster(2., notype()));
}*/

BOOST_AUTO_TEST_SUITE_END();
