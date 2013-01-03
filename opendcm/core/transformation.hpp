/*
    openDCM, dimensional constraint manager
    Copyright (C) 2013  Stefan Troeger <stefantroeger@gmx.net>

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


#ifndef DCM_TRANSFORMATION_H
#define DCM_TRANSFORMATION_H

#include <cmath>

namespace dcm {
namespace detail {

template<typename Kernel, int Dimension = 3>
class Transformation {

    typedef typename Kernel::number_type Scalar;
    typedef typename mpl::if_< boost::is_same<mpl::int_<Dimension>, mpl::int_<3> >,
            typename Kernel::Quaternion,
            Scalar>::type   			Rotation;
    typedef Eigen::Matrix<Scalar, Dimension, 1> 	Vector;

    Rotation   	m_rotation;
    Vector 	m_translation;
    Scalar  	m_scale;

#ifdef USE_LOGGING
    src::logger log;
#endif

private:
    template<int Dim>
    static inline typename boost::enable_if_c< Dim == 3, typename Kernel::Quaternion>::type initRotation() {
        return Kernel::Quaternion::Identity();
    };
    template<int Dim>
    static inline typename boost::enable_if_c< Dim == 2, Scalar>::type initRotation() {
        return 0.;
    };
    template<int Dim>
    static inline typename boost::enable_if_c< Dim == 3, typename Kernel::Quaternion>::type
    combineRotation(Rotation& r1, Rotation& r2) {
        return r2*r1;
    };
    template<int Dim>
    static inline typename boost::enable_if_c< Dim == 2, Scalar>::type
    combineRotation(Rotation& r1, Rotation& r2) {
        return r1+r2;
    };
    template<int Dim>
    inline typename boost::enable_if_c< Dim == 3, typename Kernel::Quaternion&>::type
    combineRotation(Rotation& r) {
        m_rotation = r*m_rotation;
        return m_rotation;
    };
    template<int Dim>
    inline typename boost::enable_if_c< Dim == 2, Scalar&>::type
    combineRotation(Rotation& r) {
        m_rotation += r;
        return m_rotation;
    };
    template<int Dim>
    inline typename boost::enable_if_c< Dim == 3, Vector&>::type rotateSingle(Vector& v) {
        v = m_rotation._transformVector(v);
        return v;
    };
    template<int Dim>
    inline typename boost::enable_if_c< Dim == 2, Vector&>::type rotateSingle(Vector& v) {
        Eigen::Matrix<Scalar,2,2> rot(std::cos(m_rotation), std::sin(m_rotation),
                                      -std::sin(m_rotation), std::cos(m_rotation));
        v = rot*v;
        return v;
    };
    template<int Dim>
    static inline typename boost::enable_if_c< Dim == 3, Rotation&>::type rotateInvert(Rotation& r) {
        r = r.conjugate();
        return r;
    };
    template<int Dim>
    static inline typename boost::enable_if_c< Dim == 2, Rotation&>::type rotateInvert(Rotation& r) {
        r *= -1;
        return r;
    };


public:

    Transformation(Rotation q = initRotation<Dimension>(),
                   Vector v = Vector::Zero(),
                   Scalar s = Scalar(1.))
        : m_rotation(q), m_translation(v), m_scale(s) {

#ifdef USE_LOGGING
        log.add_attribute("Tag", attrs::constant< std::string >("Transformation"));
#endif
    };

    //Accessing the values
    Rotation& rotation() {
        return m_rotation;
    };
    Rotation& rotation(Rotation r) {
        m_rotation = r;
        return m_rotation;
    };
    Vector& translation() {
        return m_translation;
    };
    Vector& translation(Vector v) {
        m_translation = v;
        return m_translation;
    };
    Scalar& scale() {
        return m_scale;
    };
    Scalar& scale(Scalar s) {
        m_scale = s;
        return m_scale;
    };

    //Transforming stuff
    void transform(Vector& v) {

        rotateSingle<Dimension>(v) += m_translation;
        v *= m_scale;
    };

    Vector& rotate(Vector& v) {
        rotateSingle<Dimension>(v);
        return v;
    };

    Vector& translate(Vector& v) {
        v += m_translation;
        return v;
    };

    Vector& scale(Vector& v) {
        v *= m_scale;
        return v;
    };

    Transformation<Kernel, Dimension>& invert() {
        rotateInvert<Dimension>(m_rotation);
        m_translation = rotate(m_translation) * (-m_scale);
        m_scale = 1./m_scale;
        return *this;
    };

    Transformation<Kernel, Dimension> inverse() {
        Transformation<Kernel, Dimension> trans(*this);
	trans.invert();
        return trans;
    };

    //successive transformations
    Transformation<Kernel, Dimension> operator+(Transformation<Kernel, Dimension>& sec) {
	Transformation<Kernel, Dimension> tmp(*this);
        return tmp+= sec;
    };
    Transformation<Kernel, Dimension>& operator+=(Transformation<Kernel, Dimension>& sec) {
        combineRotation<Dimension>(sec.rotation());
        sec.rotate(m_translation);
        m_translation += sec.translation()/m_scale;
        m_scale *= sec.scale();
        return *this;
    };
    Transformation<Kernel, Dimension> operator-(Transformation<Kernel, Dimension>& sec) {
	Transformation<Kernel, Dimension> tmp1(*this);
	Transformation<Kernel, Dimension> tmp2(sec);
	tmp2.invert();
        return tmp1+=tmp2;
    };

    Transformation<Kernel, Dimension>& operator-=(Transformation<Kernel, Dimension>& sec) {
        Transformation<Kernel, Dimension> tmp(sec);
	tmp.invert();
        return operator+=(tmp);
    };
};

}//detail
}//DCM

#endif //DCM_TRANSFORMATION
