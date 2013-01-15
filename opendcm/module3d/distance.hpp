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

#ifndef GCM_DISTANCE3D_H
#define GCM_DISTANCE3D_H

#include "geometry.hpp"
#include <opendcm/core/constraint.hpp>

namespace dcm {

template<typename Kernel>
struct Distance::type< Kernel, tag::point3D, tag::point3D > {

    typedef typename Kernel::number_type Scalar;
    typedef typename Kernel::VectorMap   Vector;

    Scalar value, sc_value;

    //template definition
    void setScale(Scalar scale) {
        sc_value = value*scale;
    };
    Scalar calculate(Vector& param1,  Vector& param2) {
        return (param1-param2).norm() - sc_value;
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

template<typename Kernel>
struct Distance::type< Kernel, tag::point3D, tag::plane3D > {
    typedef typename Kernel::number_type Scalar;
    typedef typename Kernel::VectorMap   Vector;

    Scalar value, sc_value, d_asin;
    Scalar dist, l;
    typename Kernel::Vector3 p2;

    //template definition
    void setScale(Scalar scale) {
        sc_value = value*scale;
    };
    Scalar calculate(Vector& param1,  Vector& param2) {
        //(p1-p2)°n / |n| - distance
        p2 = param2.head(3) + sc_value*param2.segment(3,3);
        dist = (param1.head(3)-p2).dot(param2.segment(3,3)) / param2.segment(3,3).norm();
        l = (param1.head(3) - p2).norm();
        d_asin = std::sqrt(1./(1-std::pow(dist/l,2)));
        return  std::asin(dist/l);
    };

    Scalar calculateGradientFirst(Vector& param1, Vector& param2, Vector& dparam1) {

        const Scalar d_dist = (dparam1.head(3)).dot(param2.segment(3,3)) / param2.segment(3,3).norm();
        Scalar d_alpha = d_dist*l;
        d_alpha -= dist*((param1.head(3) - p2).dot(dparam1.head(3)))/l;
        d_alpha /= std::pow(l,2);
        return d_asin*d_alpha;
    };

    Scalar calculateGradientSecond(Vector& param1, Vector& param2, Vector& dparam2) {

        const typename Kernel::Vector3 p1 = param1.head(3);
        const typename Kernel::Vector3 dp2 = dparam2.head(3) + sc_value*dparam2.segment(3,3);
        const typename Kernel::Vector3 n = param2.segment(3,3);
        const typename Kernel::Vector3 dn = dparam2.segment(3,3);

        const Scalar d_dist = (((-dp2).dot(n) + (p1-p2).dot(dn)) / n.norm() - (p1-p2).dot(n)* n.dot(dn)/std::pow(n.norm(),3));
        Scalar d_alpha = d_dist*l;
        d_alpha -= dist*((p1 - p2).dot(-dp2))/l;
        d_alpha /= std::pow(l,2);
        return d_asin*d_alpha;
    };

    void calculateGradientFirstComplete(Vector& param1, Vector& param2, Vector& gradient) {
        const typename Kernel::Vector3 n = param2.segment(3,3);
        gradient = d_asin*((n/n.norm())*l-dist*(param1.head(3)-p2)/l)/std::pow(l,2);
    };

    void calculateGradientSecondComplete(Vector& param1, Vector& param2, Vector& gradient) {

        const typename Kernel::Vector3 n = param2.segment(3,3);
        const typename Kernel::Vector3 p12 = (param1.head(3)-p2);
        const typename Kernel::Vector3 d_dist = ((-sc_value*n + p12)*n.norm()-(p12.dot(n)/n.norm())*n)/std::pow(n.norm(),2);
        const typename Kernel::Vector3 d_l = (-p12*sc_value)/p12.norm();

        gradient.head(3) = d_asin*((-n/n.norm())*l+dist*p12/l)/std::pow(l,2);
        gradient.segment(3,3) = d_asin * (d_dist*l - dist*d_l)/std::pow(l,2);
    };
};

template<typename Kernel>
struct Distance::type< Kernel, tag::line3D, tag::line3D > {

    typedef typename Kernel::number_type Scalar;
    typedef typename Kernel::VectorMap   Vector;
    typedef typename Kernel::Vector3     Vector3;

    Scalar value, sc_value;

    //template definition
    void setScale(Scalar scale) {
        sc_value = value*scale;
    };
    Scalar calculate(Vector& param1,  Vector& param2) {
        //diff = point1 - point2
        const Vector3 diff = param1.template head<3>() - param2.template head<3>();
        return (diff - diff.dot(param1.template segment<3>(3))*param1.template segment<3>(3)).norm() - sc_value;
    };

    Scalar calculateGradientFirst(Vector& param1, Vector& param2, Vector& dparam1) {
        const Vector3 diff = param1.template head<3>() - param2.template head<3>();
        const Vector3 dp1 = dparam1.template head<3>();
        const Vector3 p1  = param1.template head<3>();
        const Vector3 n1  = param1.template segment<3>(3);
        const Vector3 dn1 = dparam1.template segment<3>(3);
        const Vector3 p2  = param2.template head<3>();

        const Vector3 r = diff - diff.dot(n1)*n1;
        const Vector3 dr =  dp1 - (dp1.dot(n1) + diff.dot(dn1))*n1 - diff.dot(n1)*dn1;
        return r.dot(dr)/r.norm();

    };

    Scalar calculateGradientSecond(Vector& param1, Vector& param2, Vector& dparam2) {
        const Vector3 diff = param1.template head<3>() - param2.template head<3>();
        const Vector3 dp2 = dparam2.template head<3>();
        const Vector3 n1  = param1.template segment<3>(3);
        const Vector3 r = diff - diff.dot(n1)*n1;
        const Vector3 dr = -dp2 + (dp2.dot(n1))*n1;
        return r.dot(dr)/r.norm();
    };

    void calculateGradientFirstComplete(Vector& param1, Vector& param2, Vector& gradient) {
        const Vector3 diff = param1.template head<3>() - param2.template head<3>();
        const Vector3 n1  = param1.template segment<3>(3);
        const Vector3 r = diff - diff.dot(param1.template segment<3>(3))*param1.template segment<3>(3);
        gradient.template head<3>() = 2*(r - param1.template segment<3>(3)*(param1.template segment<3>(3).dot(r)));
        gradient.template segment<3>(3) = 2*((-diff)* n1.dot(r) - (diff.dot(n1))*r);
        gradient(6) = 0; //radius has nothin to do with coincidents

    };

    void calculateGradientSecondComplete(Vector& param1, Vector& param2, Vector& gradient) {
        const Vector3 diff = param1.template head<3>() - param2.template head<3>();
        const Vector3 n1  = param1.template segment<3>(3);
        const Vector3 r = diff - diff.dot(param1.template segment<3>(3))*param1.template segment<3>(3);
        gradient.template head<3>() = 2*(-r + param1.template segment<3>(3)*(param1.template segment<3>(3).dot(r)));
        gradient.template tail<4>().setZero();
    };
};

template<typename Kernel>
struct Distance::type< Kernel, tag::cylinder3D, tag::cylinder3D > : public Distance::type< Kernel, tag::line3D, tag::line3D > {};

}

#endif //GCM_DISTANCE3D_H
