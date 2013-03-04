/*
    openDCM, dimensional raint manager
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
#include <Eigen/Core>
#include <Eigen/Geometry>

#include <boost/mpl/if.hpp>

namespace dcm {
namespace detail {

template<typename Scalar, int Dim>
class Transform {

public:
    typedef Eigen::Matrix<Scalar, Dim, 1> Vector;
    typedef typename boost::mpl::if_c< Dim == 3,
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
    Transform() : m_rotation(Rotation::Identity()), m_translation(Translation::Identity()), m_scale(Scalar(1.)) {
        m_rotation.normalize();
    };
	Transform(Rotation& q) : m_rotation(q), m_translation(Translation::Identity()), m_scale(Scalar(1.)) {
        m_rotation.normalize();
    };
	Transform(Rotation& q, Translation& t) : m_rotation(q), m_translation(t), m_scale(Scalar(1.)) {
        m_rotation.normalize();
    };
	Transform(Rotation& q, Translation& t, Scalar& s) : m_rotation(q), m_translation(t), m_scale(s) {
        m_rotation.normalize();
    };
	Transform(Rotation& q, Translation& t, Scaling& s) : m_rotation(q), m_translation(t), m_scale(s.factor()) {
        m_rotation.normalize();
    };

    //access the single parts and manipulate them
    //***********************
	Rotation& rotation() {
        return m_rotation;
    }
    template<typename Derived>
    Transform& rotate( Eigen::RotationBase<Derived,Dim>& rotation) {
        m_rotation = rotation.derived().normalized()*m_rotation;
        return *this;
    }

	Translation& translation() {
        return m_translation;
    }
    Transform& translate( Translation& translation) {
        m_translation = m_translation*translation;
        return *this;
    }

	Scalar& scaling() {
        return m_scale;
    }
    Transform& scale( Scalar& scaling) {
        m_scale *= scaling;
        return *this;
    }
    Transform& scale( Scaling& scaling) {
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

    inline Transform& operator=( Translation& t) {
        m_translation = t;
        m_rotation = Rotation::Identity();
        m_scale = 1;
        return *this;
    }
    inline Transform operator*( Translation& t)  {
        Transform res = *this;
        res.translate(t);
        return res;
    }
    inline Transform& operator*=( Translation& t) {
        return translate(t);
    }

    inline Transform& operator=( Scaling& s) {
        m_scale = s.factor();
        m_translation = Translation::Identity();
        m_rotation = Rotation::Identity();
        return *this;
    }
    inline Transform operator*( Scaling& s)  {
        Transform res = *this;
        res.scale(s.factor());
        return res;
    }
    inline Transform& operator*=( Scaling& s) {
        return scale(s.factor());
    }

    template<typename Derived>
    inline Transform& operator=( Eigen::RotationBase<Derived,Dim>& r) {
        m_rotation = r.derived();
        m_rotation.normalize();
        m_translation = Translation::Identity();
        m_scale = 1;
        return *this;
    }
    template<typename Derived>
    inline Transform operator*( Eigen::RotationBase<Derived,Dim>& r)  {
        Transform res = *this;
        res.rotate(r.derived());
        return res;
    }
    template<typename Derived>
    inline Transform& operator*=( Eigen::RotationBase<Derived,Dim>& r) {
        return rotate(r.derived());
    }

    inline Transform operator* ( Transform& other)   {
        Transform res(*this);
        res*= other;
        return res;
    }
    inline Transform& operator*= ( Transform& other) {
        rotate(other.rotation());
        other.rotate(m_translation.vector());
        m_translation.vector() += other.translation().vector()/m_scale;
        m_scale *= other.scaling();
        return *this;
    }

    //transform Vectors
    //*****************
    template<typename Derived>
    inline Derived& rotate(Eigen::MatrixBase<Derived>& vec)  {
        vec = m_rotation*vec;
        return vec.derived();
    }
    template<typename Derived>
    inline Derived& translate(Eigen::MatrixBase<Derived>& vec)  {
        vec = m_translation*vec;
        return vec.derived();
    }
    template<typename Derived>
    inline Derived& scale(Eigen::MatrixBase<Derived>& vec)  {
        vec*=m_scale;
        return vec.derived();
    }
    template<typename Derived>
    inline Derived& transform(Eigen::MatrixBase<Derived>& vec)  {
        vec = (m_rotation*vec + m_translation.vector())*m_scale;
        return vec.derived();
    }
    template<typename Derived>
    inline Derived operator*( Eigen::MatrixBase<Derived>& vec)  {
        return (m_rotation*vec + m_translation.vector())*m_scale;
    }
    template<typename Derived>
    inline void operator()(Eigen::MatrixBase<Derived>& vec)  {
        transform(vec);
    }

    //Stuff
    //*****
    bool isApprox(Transform& other, Scalar prec) {
        return m_rotation.isApprox(other.rotation(), prec)
               && ((m_translation.vector()- other.translation().vector()).norm() < prec)
               && (std::abs(m_scale-other.scaling()) < prec);
    };
    void setIdentity() {
        m_rotation.setIdentity();
        m_translation = Translation::Identity();
        m_scale = 1.;
    }
    static Transform Identity() {
        return Transform();
    }

    Transform& normalize() {
        m_rotation.normalize();
        return *this;
    }
};

template<typename Scalar, int Dim>
class DiffTransform : public Transform<Scalar, Dim> {

    typedef typename Transform<Scalar, Dim>::Rotation Rotation;
    typedef typename Transform<Scalar, Dim>::Translation Translation;
    typedef Eigen::Matrix<Scalar, Dim, 3*Dim> DiffMatrix;
    
    DiffMatrix m_diffMatrix;

public:
	DiffTransform() : Transform<Scalar, Dim>() {
		m_diffMatrix.setZero();
    };
	DiffTransform(Rotation& q) : Transform<Scalar, Dim>(q) {
        m_diffMatrix.setZero();
    };
	DiffTransform(Rotation& q, Translation& t) : Transform<Scalar, Dim>(q,t) {
        m_diffMatrix.setZero();
    };
	DiffTransform(Rotation& q, Translation& t, Scalar& s) : Transform<Scalar, Dim>(q,t,s) {
        m_diffMatrix.setZero();
    };
	DiffTransform(Rotation& q, Translation& t, Scaling& s) : Transform<Scalar, Dim>(q,t,s) {
        m_diffMatrix.setZero();
    };

    DiffTransform(Transform<Scalar, Dim>& trans)
        : Transform<Scalar, Dim>(trans.rotation(), trans.translation(), trans.scaling()) {

        m_diffMatrix.setZero();
    };
    
     DiffMatrix& differential() {return m_diffMatrix;};
    inline Scalar& operator()(int f, int s) {return m_diffMatrix(f,s);};
    inline Scalar& at(int f, int s) {return m_diffMatrix(f,s);};
};

}//detail
}//DCM

/*When you overload a binary operator as a member function of a class the overload is used
 * when the first operand is of the class type.For stream operators, the first operand
 * is the stream and not (usually) the custom class.
*/
template<typename Kernel, int Dim>
std::ostream& operator<<(std::ostream& os,  dcm::detail::Transform<Kernel, Dim>& t) {
    os << "Rotation:    " << t.rotation().coeffs().transpose() << std::endl
       << "Translation: " << t.translation().vector().transpose() <<std::endl
       << "Scale:       " << t.scaling();
    return os;
}

template<typename Kernel, int Dim>
std::ostream& operator<<(std::ostream& os, dcm::detail::DiffTransform<Kernel, Dim>& t) {
    os << "Rotation:    " << t.rotation().coeffs().transpose() << std::endl
       << "Translation: " << t.translation().vector().transpose() <<std::endl
       << "Scale:       " << t.scaling() << std::endl
       << "Differential:" << std::endl<<t.differential();
    return os;
}


#endif //DCM_TRANSFORMATION
