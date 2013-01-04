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
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Geometry>

namespace dcm {
namespace detail {

template<typename Kernel, int Dim>
class Transform {

public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW_IF_VECTORIZABLE_FIXED_SIZE(typename Kernel::number_type,(Dim+1)*(Dim+1))

    typedef typename Kernel::number_type  Scalar;
    typedef Eigen::Matrix<Scalar, Dim, 1> Vector;
    typedef typename mpl::if_c< Dim== 3,
            Eigen::Quaternion<Scalar>,
            Eigen::Rotation2D<Scalar> >::type     Rotation;
    typedef Eigen::Translation<Scalar, Dim>	  Translation;
    typedef Eigen::UniformScaling<Scalar> 	  Scaling;
    typedef typename Rotation::RotationMatrixType RotationMatrix;

protected:
    Rotation   	m_rotation;
    Translation	m_translation;
    Scalar  	m_scale;

public:
    Transform(Rotation q = Rotation::Identity(),
              Translation v = Translation::Identity(),
              Scalar s = Scalar(1.))
        : m_rotation(q), m_translation(v), m_scale(s) {
	  
	  m_rotation.normalize();
    };

    //access the single parts and manipulate them
    //***********************
    const Rotation& rotation() const {
        return m_rotation;
    }
    template<typename RotationType>
    Transform& rotate(const RotationType& rotation) {
        m_rotation = rotation.normalized()*m_rotation;
        return *this;
    }

    const Translation& translation() const {
        return m_translation;
    }
    Transform& translate(const Translation& translation) {
        m_translation = m_translation*translation;
        return *this;
    }

    const Scalar& scaling() const {
        return m_scale;
    }
    Transform& scale(const Scalar& scaling) {
        m_scale *= scaling;
        return *this;
    }
    Transform& scale(const Scaling& scaling) {
        m_scale *= scaling.factor();
        return *this;
    }

    Transform& invert() {
        m_rotation = m_rotation.inverse();
        m_translation.vector() = (m_rotation*m_translation.vector()) * (-m_scale);
        m_scale = 1./m_scale;
        return *this;
    };
    Transform inverse() {
        Transform res(*this);
        res.invert();
        return res;
    };

    //operators for value manipulation
    //********************************

    inline Transform& operator=(const Translation& t) {
        m_translation = t;
        m_rotation = Rotation::Identity();
        m_scale = 1;
        return *this;
    }
    inline Transform operator*(const Translation& t) const {
        Transform res = *this;
        res.translate(t.vector());
        return res;
    }
    inline Transform& operator*=(const Translation& t) {
        return translate(t);
    }

    inline Transform& operator=(const Scaling& s) {
        m_scale = s.factor();
        m_translation = Translation::Identity();
        m_rotation = Rotation::Identity();
        return *this;
    }
    inline Transform operator*(const Scaling& s) const {
        Transform res = *this;
        res.scale(s.factor());
        return res;
    }
    inline Transform& operator*=(const Scaling& s) {
        return scale(s.factor());
    }

    template<typename Derived>
    inline Transform& operator=(const Eigen::RotationBase<Derived,Dim>& r) {
        m_rotation = r.derived();
	m_rotation.normalize();
        m_translation = Translation::Identity();
        m_scale = 1;
        return *this;
    }
    template<typename Derived>
    inline Transform operator*(const Eigen::RotationBase<Derived,Dim>& r) const {
        Transform res = *this;
        res.rotate(r.derived());
        return res;
    }
    template<typename Derived>
    inline Transform& operator*=(const Eigen::RotationBase<Derived,Dim>& r) {
        return rotate(r.derived());
    }

    inline Transform operator* (const Transform& other) const  {
        Transform res(*this);
        res*= other;
        return res;
    }
    inline Transform& operator*= (const Transform& other) {
        rotate(other.rotation());
        m_translation.vector() = other.rotate(m_translation.vector());
        m_translation.vector() += other.translation().vector()/m_scale;
        m_scale *= other.scaling();
        return *this;
    }

    //transform Vectors
    //*****************
    inline Vector rotate(const Vector& vec) const {
        return m_rotation*vec;
    }
    inline Vector translate(const Vector& vec) const {
        return m_translation*vec;
    }
    inline Vector scale(const Vector& vec) const {
        return vec*m_scale;
    }
    inline Vector& transform(Vector& vec) const {
        vec = (m_rotation*vec + m_translation.vector())*m_scale;
        return vec;
    }
    inline Vector operator*(const Vector& vec) const {
        return (*this)(vec);
    }
    inline Vector operator()(const Vector& vec) const {
        return (m_rotation*vec + m_translation.vector())*m_scale;
    }

    //Stuff
    //*****
    bool isApprox(const Transform& other, Scalar prec) const {
        return m_rotation.isApprox(other.rotation(), prec)
               && m_translation.isApprox(other.translation(), prec)
               && (std::abs(m_scale-other.scaling()) < prec);
    };
    void setIdentity() {
        m_rotation.setIdentity();
        m_translation = Translation::Identity();
        m_scale = 1.;
    }
    static const Transform Identity() {
        return Transform(Rotation::Identity(), Translation::Identity, 1.);
    }


    /*

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
    };*/
};

}//detail
}//DCM

#endif //DCM_TRANSFORMATION
