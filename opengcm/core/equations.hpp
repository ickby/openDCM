/*
    openGCM, geometric constraint manager
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

#ifndef GCM_EQUATIONS_H
#define GCM_EQUATIONS_H

#include <assert.h>

namespace gcm {

struct no_option {};

struct Distance {

    typedef double option_type;
    option_type value;

    Distance() : value(0) {};

    Distance& operator=(const option_type& val) {
        value = val;
    };

    template< typename Kernel, typename Tag1, typename Tag2 >
    struct type {

        typedef typename Kernel::number_type Scalar;
        typedef typename Kernel::VectorMap   Vector;

        Scalar value;
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
};

//the possible directions
enum Direction { Same, Opposite, Both };

struct Parallel {

    typedef Direction option_type;
    option_type value;

    Parallel() : value(Both) {};

    Parallel& operator=(const option_type& val) {
        value = val;
    };

    template< typename Kernel, typename Tag1, typename Tag2 >
    struct type {

        typedef typename Kernel::number_type Scalar;
        typedef typename Kernel::VectorMap   Vector;
	
	option_type value;

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
};

struct Angle {

    typedef double option_type;
    option_type value;

    Angle() : value(0) {};

    Parallel& operator=(const option_type& val) {
        value = val;
    };

    template< typename Kernel, typename Tag1, typename Tag2 >
    struct type {

        typedef typename Kernel::number_type Scalar;
        typedef typename Kernel::VectorMap   Vector;

	option_type value;
	
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
};

//static is needed to restrain the scope of the objects to the current compilation unit. Without it
//every compiled file including this header would define these as global and the linker would find
//multiple definitions of the same objects
static Distance distance;
static Parallel parallel;
static Angle    angle;

};

#endif //GCM_EQUATIONS_H

