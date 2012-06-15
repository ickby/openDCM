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
    GNU General Public License for more detemplate tails.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/

#ifndef DCM_CONSTRAINT3D_H
#define DCM_CONSTRAINT3D_H

#include "geometry.hpp"

namespace dcm {


template<typename Kernel, typename Tag1, typename Tag2>
struct Coincident3D {

    typedef typename Kernel::number_type Scalar;
    typedef typename Kernel::VectorMap   Vector;

    Scalar calculate(Vector& param1,  Vector& param2)  {
        assert(false);
    };

    Scalar calculateGradientFirst(Vector& param1,  Vector& param2, Vector& dparam1) {
        assert(false);
    };

    Scalar calculateGradientSecond(Vector& param1,  Vector& param2, Vector& dparam2)  {
        assert(false);
    };

    void calculateGradientFirstComplete(Vector& param1,  Vector& param2, Vector& gradient) {
        assert(false);
    };

    void calculateGradientSecondComplete(Vector& param1,  Vector& param2, Vector& gradient) {
        assert(false);
    };

};

/*******************************************************************************
 * 			Distance Constraint
 * *****************************************************************************
 */

template< typename Kernel, typename Tag1, typename Tag2 >
struct Distance3D {

    typedef typename Kernel::number_type Scalar;
    typedef typename Kernel::VectorMap   Vector;
    Scalar m_distance;

    Distance3D(Scalar d = 0) : m_distance(d) {};

    //template definition
    Scalar calculate(Vector& param1,  Vector& param2) {
        assert(false);
    };
    Scalar calculateGradientFirst(Vector& param1, Vector& param2, Vector& dparam1) {
        assert(false);
    };
    Scalar calculateGradientSecond(Vector& param1, Vector& param2, Vector& dparam2) {
        assert(false);
    };
    void calculateGradientFirstComplete(Vector& param1, Vector& param2, Vector& gradient) {
        assert(false);
    };
    void calculateGradientSecondComplete(Vector& param1, Vector& param2, Vector& gradient) {
        assert(false);
    };
};

template< typename Kernel >
struct Distance3D< Kernel, tag::point3D, tag::point3D > {

    typedef typename Kernel::number_type Scalar;
    typedef typename Kernel::VectorMap   Vector;

    Scalar m_distance;

    Distance3D(Scalar d = 0) : m_distance(d) {};

    //template definition
    Scalar calculate(Vector& param1,  Vector& param2) {
        return (param1-param2).norm() - m_distance;
    };

    Scalar calculateGradientFirst(Vector& param1, Vector& param2, Vector& dparam1) {
        return (param1-param2).dot(dparam1) / (param1-param2).norm();
    };

    Scalar calculateGradientSecond(Vector& param1, Vector& param2, Vector& dparam2) {
        return (param1-param2).dot(-dparam2) / (param1-param2).norm();
    };

    void calculateGradientFirstComplete(Vector& param1, Vector& param2, Vector& gradient) {
        gradient = (param1-param2) / (param1-param2).norm();
    };

    void calculateGradientSecondComplete(Vector& param1, Vector& param2, Vector& gradient) {
        gradient = (param2-param1) / (param1-param2).norm();
    };
};

//remember: only valid for parallel planes (as intersecting have always minimal distance 0) 
template< typename Kernel >
struct Distance3D< Kernel, tag::plane3D, tag::plane3D > {

    typedef typename Kernel::number_type Scalar;
    typedef typename Kernel::VectorMap   Vector;

    Scalar m_distance;

    Distance3D(Scalar d = 0) : m_distance(d) {};

    //template definition
    Scalar calculate(Vector& param1,  Vector& param2) {
	//(p1-p2)°n / |n| - distance
	return  (param1.head(3)-param2.head(3)).dot(param2.tail(3)) / param2.tail(3).norm() - m_distance;
    };

    Scalar calculateGradientFirst(Vector& param1, Vector& param2, Vector& dparam1) {
      //dp1°n / |n|
        return (dparam1.head(3)).dot(param2.tail(3)) / param2.tail(3).norm();
    };

    Scalar calculateGradientSecond(Vector& param1, Vector& param2, Vector& dparam2) {
	const typename Kernel::Vector3 p1 = param1.head(3);
	const typename Kernel::Vector3 p2 = param2.head(3);
	const typename Kernel::Vector3 dp2 = dparam2.head(3);
	const typename Kernel::Vector3 n = param2.tail(3);
	const typename Kernel::Vector3 dn = dparam2.tail(3);
	
        return ((-dp2).dot(n) + (p1-p2).dot(dn)) / n.norm() - (p1-p2).dot(n)* n.dot(dn)/std::pow(n.norm(),3);
    };

    void calculateGradientFirstComplete(Vector& param1, Vector& param2, Vector& gradient) {
        gradient = param2.tail(3) / param2.tail(3).norm();
    };

    void calculateGradientSecondComplete(Vector& param1, Vector& param2, Vector& gradient) {
	const typename Kernel::Vector3 p1m2 = param1.head(3) - param2.head(3);
	const typename Kernel::Vector3 n = param2.tail(3);
	
        gradient.head(3) = -n / n.norm();
	gradient.tail(3) = (p1m2)/n.norm() - (p1m2).dot(n)*n/std::pow(n.norm(),3);
    };
};


/*******************************************************************************
 * 			Parallel Constraint
 * *****************************************************************************
 */

//the possible directions
enum Direction { Same, Opposite };

//the calculations( same as we always calculate directions we can outsource the work to this functions)
namespace parallel {

template<typename Kernel, typename T>
inline typename Kernel::number_type calc(T d1,
        T d2,
        Direction dir)  {

    switch(dir) {
        case Same:
            return (d1-d2).norm();
        case Opposite:
            return (d1+d2).norm();
    }
};


template<typename Kernel, typename T>
inline typename Kernel::number_type calcGradFirst(T d1,
        T d2,
        T dd1,
        Direction dir)  {

    switch(dir) {
        case Same:
            return (d1-d2).dot(dd1) / (d1-d2).norm();
        case Opposite:
            return (d1+d2).dot(dd1) / (d1+d2).norm();
    }
};

template<typename Kernel, typename T>
inline typename Kernel::number_type calcGradSecond(T d1,
        T d2,
        T dd2,
        Direction dir)  {

    switch(dir) {
        case Same:
            return (d1-d2).dot(-dd2) / (d1-d2).norm();
        case Opposite:
            return (d1+d2).dot(dd2) / (d1+d2).norm();
    }
};

template<typename Kernel, typename T>
inline void calcGradFirstComp(T d1,
                              T d2,
                              T grad,
                              Direction dir)  {

    switch(dir) {
        case Same:
            grad = (d1-d2) / (d1-d2).norm();
            return;
        case Opposite:
            grad = (d1+d2) / (d1+d2).norm();
            return;
    }
};

template<typename Kernel, typename T>
inline void calcGradSecondComp(T d1,
                               T d2,
                               T grad,
                               Direction dir)  {

    switch(dir) {
        case Same:
            grad = (d2-d1) / (d1-d2).norm();
            return;
        case Opposite:
            grad = (d2+d1) / (d1+d2).norm();
            return;
    }
};

}

template< typename Kernel, typename Tag1, typename Tag2 >
struct Parallel3D {

    typedef typename Kernel::number_type Scalar;
    typedef typename Kernel::VectorMap   Vector;
    Direction m_dir;

    Parallel3D(Direction d = Same) : m_dir(d) {
      Base::Console().Message("choosen direction (0=same, 1=opposite): %d\n",m_dir);
    };

    //template definition
    Scalar calculate(Vector& param1,  Vector& param2) {
        assert(false);
    };
    Scalar calculateGradientFirst(Vector& param1, Vector& param2, Vector& dparam1) {
        assert(false);
    };
    Scalar calculateGradientSecond(Vector& param1, Vector& param2, Vector& dparam2) {
        assert(false);
    };
    void calculateGradientFirstComplete(Vector& param1, Vector& param2, Vector& gradient) {
        assert(false);
    };
    void calculateGradientSecondComplete(Vector& param1, Vector& param2, Vector& gradient) {
        assert(false);
    };
};

template< typename Kernel >
struct Parallel3D< Kernel, tag::line3D, tag::line3D > {

    typedef typename Kernel::number_type Scalar;
    typedef typename Kernel::VectorMap   Vector;

    Direction m_dir;

    Parallel3D(Direction d = Same) : m_dir(d) {};

    //template definition
    Scalar calculate(Vector& param1,  Vector& param2) {
        return parallel::calc<Kernel>(param1.template tail<3>(), param2.template tail<3>(), m_dir);
    };
    Scalar calculateGradientFirst(Vector& param1, Vector& param2, Vector& dparam1) {
        return parallel::calcGradFirst<Kernel>(param1.template tail<3>(), param2.template tail<3>(), dparam1.template tail<3>(), m_dir);
    };
    Scalar calculateGradientSecond(Vector& param1, Vector& param2, Vector& dparam2) {
        return parallel::calcGradSecond<Kernel>(param1.template tail<3>(), param2.template tail<3>(), dparam2.template tail<3>(), m_dir);
    };
    void calculateGradientFirstComplete(Vector& param1, Vector& param2, Vector& gradient) {
        gradient.template head<3>().setZero();
        parallel::calcGradFirstComp<Kernel>(param1.template tail<3>(), param2.template tail<3>(), gradient.template tail<3>(), m_dir);
    };
    void calculateGradientSecondComplete(Vector& param1, Vector& param2, Vector& gradient) {
        gradient.template head<3>().setZero();
        parallel::calcGradSecondComp<Kernel>(param1.template tail<3>(), param2.template tail<3>(), gradient.template tail<3>(), m_dir);
    };
};

//planes like lines have the direction as segment 3-5, so we can use the same implementations
template< typename Kernel >
struct Parallel3D< Kernel, tag::plane3D, tag::plane3D > : public Parallel3D<Kernel, tag::line3D, tag::line3D> {
  Parallel3D(Direction d = Same) : Parallel3D<Kernel, tag::line3D, tag::line3D>(d) {};
};
template< typename Kernel >
struct Parallel3D< Kernel, tag::line3D, tag::plane3D > : public Parallel3D<Kernel, tag::line3D, tag::line3D> {
  Parallel3D(Direction d = Same) : Parallel3D<Kernel, tag::line3D, tag::line3D>(d) {};
};


/*******************************************************************************
 * 			Angle Constraint
 * *****************************************************************************
 */

//the calculations( same as we always calculate directions we can outsource the work to this functions)
namespace angle {

template<typename Kernel, typename T>
inline typename Kernel::number_type calc(T d1,
        T d2,
        typename Kernel::number_type angle)  {

    return d1.dot(d2) / (d1.norm()*d2.norm()) - angle;
};


template<typename Kernel, typename T>
inline typename Kernel::number_type calcGradFirst(T d1,
        T d2,
        T dd1)  {

    typename Kernel::number_type norm = d1.norm()*d2.norm();
    return  dd1.dot(d2)/norm - (d1.dot(d2) * (d1.dot(dd1)/d1.norm())*d2.norm()) / std::pow(norm,2);
};

template<typename Kernel, typename T>
inline typename Kernel::number_type calcGradSecond(T d1,
        T d2,
        T dd2)  {

    typename Kernel::number_type norm = d1.norm()*d2.norm();
    return  d1.dot(dd2)/norm - (d1.dot(d2) * (d2.dot(dd2)/d2.norm())*d1.norm()) / std::pow(norm,2);
};

template<typename Kernel, typename T>
inline void calcGradFirstComp(T d1,
                              T d2,
                              T grad)  {

    typename Kernel::number_type norm = d1.norm()*d2.norm();
    grad = d2/norm - (d1.dot(d2)/std::pow(norm,2))*d1/d1.norm();
};

template<typename Kernel, typename T>
inline void calcGradSecondComp(T d1,
                               T d2,
                               T grad)  {

    typename Kernel::number_type norm = d1.norm()*d2.norm();
    grad = d1/norm - (d1.dot(d2)/std::pow(norm,2))*d2/d2.norm();
};

}

template< typename Kernel, typename Tag1, typename Tag2 >
struct Angle3D {

    typedef typename Kernel::number_type Scalar;
    typedef typename Kernel::VectorMap   Vector;
    Scalar m_angle;

    Angle3D(Scalar d = 0.) : m_angle(std::cos(d)) {};

    //template definition
    Scalar calculate(Vector& param1,  Vector& param2) {
        assert(false);
    };
    Scalar calculateGradientFirst(Vector& param1, Vector& param2, Vector& dparam1) {
        assert(false);
    };
    Scalar calculateGradientSecond(Vector& param1, Vector& param2, Vector& dparam2) {
        assert(false);
    };
    void calculateGradientFirstComplete(Vector& param1, Vector& param2, Vector& gradient) {
        assert(false);
    };
    void calculateGradientSecondComplete(Vector& param1, Vector& param2, Vector& gradient) {
        assert(false);
    };
};

template< typename Kernel >
struct Angle3D< Kernel, tag::line3D, tag::line3D > {

    typedef typename Kernel::number_type Scalar;
    typedef typename Kernel::VectorMap   Vector;

    Scalar m_angle;

    Angle3D(Scalar d = 0.) : m_angle(std::cos(d)) {};

    //template definition
    Scalar calculate(Vector& param1,  Vector& param2) {
        return angle::calc<Kernel>(param1.template tail<3>(), param2.template tail<3>(), m_angle);
    };
    Scalar calculateGradientFirst(Vector& param1, Vector& param2, Vector& dparam1) {
        return angle::calcGradFirst<Kernel>(param1.template tail<3>(), param2.template tail<3>(), dparam1.template tail<3>());
    };
    Scalar calculateGradientSecond(Vector& param1, Vector& param2, Vector& dparam2) {
        return angle::calcGradSecond<Kernel>(param1.template tail<3>(), param2.template tail<3>(), dparam2.template tail<3>());
    };
    void calculateGradientFirstComplete(Vector& param1, Vector& param2, Vector& gradient) {
        gradient.template head<3>().setZero();
        angle::calcGradFirstComp<Kernel>(param1.template tail<3>(), param2.template tail<3>(), gradient.template tail<3>());
    };
    void calculateGradientSecondComplete(Vector& param1, Vector& param2, Vector& gradient) {
        gradient.template head<3>().setZero();
        angle::calcGradSecondComp<Kernel>(param1.template tail<3>(), param2.template tail<3>(), gradient.template tail<3>());
    };
};

//planes like lines have the direction as segment 3-5, so we can use the same implementations
template< typename Kernel >
struct Angle3D< Kernel, tag::plane3D, tag::plane3D > : public Angle3D<Kernel, tag::line3D, tag::line3D> {};
template< typename Kernel >
struct Angle3D< Kernel, tag::line3D, tag::plane3D > : public Angle3D<Kernel, tag::line3D, tag::line3D> {};

}

#endif //DCM_CONSTRAINT3D_H
