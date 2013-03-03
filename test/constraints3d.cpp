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

#include "opendcm/core.hpp"
#include "opendcm/module3d.hpp"
#include "opendcm/modulepart.hpp"

#include <boost/test/unit_test.hpp>

typedef Eigen::Matrix<double, 3,1> point_t;
struct line_t : public Eigen::Matrix<double, 6,1> {};
struct plane_t : public Eigen::Matrix<double, 6,1> {};

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
    CheckSolver() {};

    template<typename Functor>
    int solve(typename Kernel::MappedEquationSystem& sys, Functor& f) {

        int jcount = 1000;
        sys.recalculate();
        double val = sys.Residual(0);

        double diffmax = 0;
        int fails = 0;
        for(int i=0; i<sys.Parameter.rows(); i++) {

            //varying parameter i
            for(int j=0; j<jcount; j++) {

                //get new value
                sys.Parameter(i) += 1e-3;
                sys.recalculate();

                //get the residual diff between the last two calculations
                double rdiff = (sys.Residual(0) - val)/1e-3;

                // if(sys.Jacobi(0,i) > 2.)
                //     std::cout<<"huge differential: "<<sys.Residual(0)<<std::endl;

                diffmax = std::max(sys.Jacobi(0,i), diffmax);

                //compare with the algeraic caclculated diff
                if(std::abs(sys.Jacobi(0,i) - rdiff) > 1e-2) {
                    std::cout<<"iter: "<<j<<", parameter: "<<i<<", calc: "<<sys.Jacobi(0,i)<<" vs real: "<<rdiff<<std::endl;
                    fails++;
                    if(fails>3)
                        throw(std::exception()); //and throw if the difference is too much
                }
                val = sys.Residual(0);
            };
        };
        return 1;
    }
};



using namespace dcm;

BOOST_AUTO_TEST_SUITE(constraint3d_test_suit);

typedef dcm::Kernel<double, CheckSolver> Kernel;
typedef Module3D< mpl::vector3<point_t, line_t, plane_t > > Module;
typedef dcm::ModulePart< mpl::vector1< place > > ModulePart;
typedef System<Kernel, Module::type, ModulePart::type> System;

typedef Module::type<System>::Geometry3D geom;
typedef boost::shared_ptr<geom> geom_ptr;

typedef Module::type<System>::Constraint3D cons;
typedef boost::shared_ptr<cons> cons_ptr;

typedef ModulePart::type<System>::Part part;
typedef boost::shared_ptr<part> part_ptr;

typedef Module::type<System>::vertex_prop vertex_prop;
typedef System::Cluster ClusterGraph;

//automated checking of differentials
template<typename T1, typename T2, typename CT>
struct constraint_checker {

    System& system;

    constraint_checker(System& sys) : system(sys) { };

    template<typename T>
    bool check_normal(T val) {

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
        type = val;
        cons_ptr c = system.createConstraint3D(g1, g2, type);

        try {
            system.solve();
            system.removeGeometry3D(g1);
            system.removeGeometry3D(g2);
            return true;
        } catch(...) {
            system.removeGeometry3D(g1);
            system.removeGeometry3D(g2);
            return false;
        }
    }

    template<typename T>
    bool check_cluster(T val) {

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
        type = val;
        cons_ptr c = system.createConstraint3D(g1, g2, type);

        try {
            system.solve();
            system.removePart(p1);
            system.removePart(p2);
            return true;
        } catch(...) {
            system.removePart(p1);
            system.removePart(p2);
            return false;
        }


    }

};


BOOST_AUTO_TEST_CASE(constraint3d_distance) {

    System sys;
    constraint_checker<point_t, point_t, dcm::Distance> checker(sys);
    BOOST_REQUIRE(checker.check_normal(2));
    BOOST_REQUIRE(checker.check_cluster(2));

    constraint_checker<point_t, plane_t, dcm::Distance> checker2(sys);
    BOOST_REQUIRE(checker2.check_normal(2));
    BOOST_REQUIRE(checker2.check_cluster(2));

    constraint_checker<point_t, line_t, dcm::Distance> checker3(sys);
    BOOST_REQUIRE(checker3.check_normal(2));
    BOOST_REQUIRE(checker3.check_cluster(2));
}


BOOST_AUTO_TEST_SUITE_END();
