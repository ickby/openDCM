/*
    openDCM, dimensional constraint manager
    Copyright (C) 2013  Stefan Troeger <stefantroeger@gmx.net>

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


#ifndef DCM_TRANSFORMATION_H
#define DCM_TRANSFORMATION_H

#include <cmath>
#include <Eigen/Core>
#include <Eigen/Geometry>

#include <boost/mpl/if.hpp>

namespace dcm {
namespace details {
    
//it is not possible to access Derived typedefs in TransformBase global scope (only function level)
//hence we provide this traits class in which all transform subclasses need to bring their definitions 
   
template<typename Transform>
struct transform_traits {};

template<typename Derived>
class TransformBase {

public:
    static const int Dimension = transform_traits<Derived>::Dimension;
    typedef typename transform_traits<Derived>::Scalar              Scalar;
    typedef Eigen::Matrix<Scalar, Dimension, 1>   Vector;
    typedef typename transform_traits<Derived>::Rotation            Rotation;
    typedef typename transform_traits<Derived>::Translation         Translation;
    typedef typename transform_traits<Derived>::RotationMatrix      RotationMatrix;

public:
    TransformBase();

    TransformBase(const Rotation& r);
    TransformBase(const Translation& t);
    TransformBase(const Rotation& r, const Translation& t);

    inline Derived& derived() {return *static_cast<Derived*>(this);};
    inline const Derived& derived() const {return *static_cast<const Derived*>(this);};
        
    //access the single parts and manipulate them
    //***********************
    const Rotation& rotation() const;
    template<typename OtherDerived>
    Derived& setRotation(const Eigen::RotationBase<OtherDerived,Dimension>& rotation);
    //rotate does not transform the translational part! For this use multiplication operators
    template<typename OtherDerived>
    Derived& rotate(const Eigen::RotationBase<OtherDerived,Dimension>& rotation);

    const Translation& translation() const;
    Derived& setTranslation(const Translation& translation);
    Derived& translate(const Translation& translation);
   
    Eigen::Transform<Scalar, Dimension, Eigen::AffineCompact> transformation();

    Derived& invert();
    Derived  inverse() const;

    //operators for value manipulation
    //********************************

    Derived& operator=(const Translation& t);
    Derived operator*(const Translation& s) const;
    Derived& operator*=(const Translation& t);

    template<typename OtherDerived>
    Derived& operator=(const Eigen::RotationBase<OtherDerived,Dimension>& r);
    template<typename OtherDerived>
    Derived operator*(const Eigen::RotationBase<OtherDerived,Dimension>& r) const;
    template<typename OtherDerived>
    Derived& operator*=(const Eigen::RotationBase<OtherDerived,Dimension>& r);

    template<typename OtherDerived>
    Derived operator* (const TransformBase<OtherDerived>& other) const;
    template<typename OtherDerived>
    Derived& operator*= (const TransformBase<OtherDerived>& other);

    //transform Vectors
    //*****************
    template<typename OtherDerived>
    OtherDerived& rotate(Eigen::MatrixBase<OtherDerived>& vec) const;
    template<typename OtherDerived>
    OtherDerived& translate(Eigen::MatrixBase<OtherDerived>& vec) const;
    template<typename OtherDerived>
    OtherDerived& transform(Eigen::MatrixBase<OtherDerived>& vec) const;
    template<typename OtherDerived>
    OtherDerived operator*(const Eigen::MatrixBase< OtherDerived >& vec) const;
    template<typename OtherDerived>
    void operator()(Eigen::MatrixBase<OtherDerived>& vec) const;

    //Stuff
    //*****
    template<typename OtherDerived>
    bool isApprox(const TransformBase<OtherDerived>& other, Scalar prec) const;
    void setIdentity();
    static const Derived Identity();

    Derived& normalize();

public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

/* When you overload a binary operator as a member function of a class the overload is used
 * when the first operand is of the class type. For stream operators, the first operand
 * is the stream and not (usually) the custom class.
*/
template<typename charT, typename traits, typename Derived>
std::basic_ostream<charT,traits>& operator<<(std::basic_ostream<charT,traits>& os, const dcm::details::TransformBase<Derived>& t);


/**
 * @brief Default Transformation class with internal storage
 * This class can be used to store and handle transformations. It has an internal storage to store the
 * transformation data. Hence it can be fully manipulated, which means it is default constructible,
 * and the multiplication operators returning Transform type work. The itnernal storage for rotations
 * is a quaternion, hence normalisation is cheap. 
 * @note As this class uses Quaternion as storage, it cannot represent a [0] rotation matrix, hence 
 *       using this for derivatives where the rotated part should become 0 is not possible.
 */
template<typename Scalar_, int Dim>
class Transform : public TransformBase<Transform<Scalar_, Dim>> {
    
public:
    typedef TransformBase<Transform<Scalar_, Dim>> Base;
    
    typedef Scalar_ Scalar;
    static const int Dimension = Dim;
    
    typedef Eigen::Quaternion<Scalar>             Rotation;
    typedef Eigen::Translation<Scalar, Dim>       Translation;
    typedef Eigen::Matrix<Scalar, 3, 1>           TranslationVector;
    typedef typename Rotation::RotationMatrixType RotationMatrix;
    
    Transform() : Base() { setIdentity();}; //we need set identity here, base can't do it as it would override the map values
    Transform(const Rotation& r) : Base(r) {};
    Transform(const Translation& t) : Base(t) {};
    Transform(const Rotation& r, const Translation& t) : Base(r,t) {};
    
    //create the functions needed from the base class
    void setIdentity() {
        m_rotation = Rotation::Identity();
        m_translation = Translation::Identity();
    };
    void normalize() {m_rotation.normalize();}
    Transform& invert() {
        m_rotation = m_rotation.inverse();
        m_translation.vector() = (m_rotation*m_translation.vector());
        return *this;
    };
    
    const Rotation& rotation() const {return m_rotation;};
    const Translation& translation() const {return m_translation;};
    const TranslationVector& translationVector() const {return m_translation.vector();};
    
    template<typename Derived>
    void setRotation(const Eigen::RotationBase<Derived, Dim>& r) {
        m_rotation = r.derived();
    };
    void setTranslation(const Translation& t) {
        m_translation = t;
    };
    
private:
    Rotation m_rotation;
    Translation m_translation;
};

template<typename S, int D>
struct transform_traits<Transform<S,D>> {
    typedef S Scalar;
    static const int Dimension = D;
    
    typedef Eigen::Quaternion<Scalar>             Rotation;
    typedef Eigen::Translation<Scalar, D>         Translation;
    typedef typename Rotation::RotationMatrixType RotationMatrix;
};


/**
 * @brief Map Transformation class without internal storage
 * This class can be used to access external transformation data and make it available as default 
 * dcm transformation. It has no own storage, but maps to external one. Hence it can save many 
 * uneeded copies if the data is already stored. Note that this type is not default constructible,
 * also as it has no internal storage the multiplication opertors don't work.
 * The rotation storage sheme of this class is a matrix.
 */
template<typename Scalar_, int Dim>
class MapMatrixTransform : public TransformBase<MapMatrixTransform<Scalar_, Dim>> {
    
public:
    typedef TransformBase<Transform<Scalar_, Dim>> Base;
    
    typedef Scalar_ Scalar;
    static const int Dimension = Dim;
    
    typedef Eigen::Map<Eigen::Matrix<Scalar, Dim, Dim>>  Rotation;
    typedef Eigen::Map<Eigen::Matrix<Scalar, Dim, 1  >>  Translation;
    typedef Eigen::Matrix<Scalar, Dim, Dim>              RotationMatrix;
       
    template<typename D1, typename D2>
    MapMatrixTransform(Eigen::MatrixBase<D1>& rot, Eigen::MatrixBase<D2>& trans) 
            : m_rotation(&rot(0,0)), m_translation(&trans(0)) {};
    
    //create the functions needed from the base class
    void setIdentity() {
        m_rotation.setIdentity();
        m_translation = Translation::Identity();
    };
    void normalize() {
        m_rotation = Eigen::Quaternion<Scalar>(m_rotation).normalize().toRotationMatrix();
    }
    
    MapMatrixTransform& invert() {
        m_rotation.transposeInPlace();
        m_translation = (m_rotation*m_translation);
        return *this;
    };
    
    const Rotation& rotation() const {return m_rotation;};
    const Translation& translation() const {return m_translation;};
    const Translation& translationVector() const {return m_translation;};
    
    template<typename Derived>
    void setRotation(const Eigen::RotationBase<Derived, Dim>& r) {
        m_rotation = r.derived().toRotationMatrix();
    };
    void setTranslation(const Translation& t) {
        m_translation = t;
    };
    
private:
    Rotation m_rotation;
    Translation m_translation;
};

template<typename S, int D>
struct transform_traits<MapMatrixTransform<S,D>> {
    typedef S Scalar;
    static const int Dimension = D;
    
    typedef Eigen::Map<Eigen::Matrix<Scalar, D, D>> Rotation;
    typedef Eigen::Map<Eigen::Matrix<Scalar, D, 1>> Translation;
    typedef Eigen::Matrix<Scalar, D, D>             RotationMatrix;
};

/**********************************************************************************************************************************
 *
 *                      IMPELEMNTATION
 * 
 * ********************************************************************************************************************************/

template<typename Derived>
TransformBase<Derived>::TransformBase() {};

template<typename Derived>
TransformBase<Derived>::TransformBase(const Rotation& r) {
    
    derived().setIdentity();
    derived().setRotation(r);
    derived().normalize();
};

template<typename Derived>
TransformBase<Derived>::TransformBase(const Translation& t) {
       
    derived().setIdentity();
    derived().setTranslation(t);
};

template<typename Derived>
TransformBase<Derived>::TransformBase(const Rotation& r, const Translation& t) {
    
    derived().setRotation(r);
    derived().setTranslation(t);
    derived().normalize();
};

template<typename Derived>
const typename TransformBase<Derived>::Rotation& TransformBase<Derived>::rotation() const {
    return derived().rotation();
}

template<typename Derived>
template<typename OtherDerived>
Derived& TransformBase<Derived>::setRotation(const Eigen::RotationBase<OtherDerived,Dimension>& rotation) {
    derived().setRotation(rotation);
    return derived();
}

template<typename Derived>
template<typename OtherDerived>
Derived& TransformBase<Derived>::rotate(const Eigen::RotationBase<OtherDerived,Dimension>& rotation) {
    derived().setRotation(rotation.derived().normalized()*derived().rotation());
    return derived();
}

template<typename Derived>
const typename TransformBase<Derived>::Translation& TransformBase<Derived>::translation() const {
    return derived().translation();
}

template<typename Derived>
Derived& TransformBase<Derived>::setTranslation(const Translation& translation) {
    derived().setTranslation(translation);
    return derived();
}

template<typename Derived>
Derived& TransformBase<Derived>::translate(const Translation& translation) {
    derived().setTranslation(derived().translation()*translation);
    return derived();
}

template<typename Derived>
Derived& TransformBase<Derived>::invert() {
    derived().invert();
    return derived();
};
template<typename Derived>
Derived TransformBase<Derived>::inverse() const {
    Derived res(derived());
    res.invert();
    return res;
};

template<typename Derived>
inline Derived& TransformBase<Derived>::operator=(const Translation& t) {
    derived().setIdentity();
    derived().translation() = t;
    return derived();
}
template<typename Derived>
inline Derived TransformBase<Derived>::operator*(const Translation& t) const {
    Derived res = derived();
    res.translate(t);
    return res;
}
template<typename Derived>
inline Derived& TransformBase<Derived>::operator*=(const Translation& t) {
    return translate(t);
}

template<typename Derived>
template<typename OtherDerived>
inline Derived& TransformBase<Derived>::operator=(const Eigen::RotationBase<OtherDerived,Dimension>& r) {
    derived().setIdentity();
    derived().setRotation(r.derived());
    derived().normalize();
    return derived();
}
template<typename Derived>
template<typename OtherDerived>
inline Derived TransformBase<Derived>::operator*(const Eigen::RotationBase<OtherDerived,Dimension>& r) const {
    Derived res = derived();
    res.rotate(r.derived());
    return res;
}
template<typename Derived>
template<typename OtherDerived>
inline Derived& TransformBase<Derived>::operator*=(const Eigen::RotationBase<OtherDerived,Dimension>& r) {
    rotate(r);
    derived().setTranslation(r*derived().translationVector());
    return derived();
}

template<typename Derived>
template<typename OtherDerived>
inline Derived TransformBase<Derived>::operator* (const TransformBase<OtherDerived>& other) const  {
    Derived res(derived());
    res*= other;
    return res;
}
template<typename Derived>
template<typename OtherDerived>
inline Derived& TransformBase<Derived>::operator*= (const TransformBase<OtherDerived>& other) {
    rotate(other.rotation());
    derived().setTranslation(Translation(other.rotation()*derived().translationVector() + other.derived().translationVector()));
    return derived();
}

template<typename Derived>
template<typename OtherDerived>
inline OtherDerived& TransformBase<Derived>::rotate(Eigen::MatrixBase<OtherDerived>& vec) const {
    vec = rotation()*vec;
    return vec.derived();
}

template<typename Derived>
template<typename OtherDerived>
inline OtherDerived& TransformBase<Derived>::translate(Eigen::MatrixBase<OtherDerived>& vec) const {
    vec = translation()*vec;
    return vec.derived();
}

template<typename Derived>
template<typename OtherDerived>
inline OtherDerived& TransformBase<Derived>::transform(Eigen::MatrixBase<OtherDerived>& vec) const {
    vec = (rotation()*vec + derived().translationVector());
    return vec.derived();
}

template<typename Derived>
template<typename OtherDerived>
inline OtherDerived TransformBase<Derived>::operator*(const Eigen::MatrixBase<OtherDerived>& vec) const {
    return (rotation()*vec + derived().translationVector());
}

template<typename Derived>
template<typename OtherDerived>
inline void TransformBase<Derived>::operator()(Eigen::MatrixBase<OtherDerived>& vec) const {
    transform(vec);
}

template<typename Derived>
template<typename OtherDerived>
bool TransformBase<Derived>::isApprox(const TransformBase<OtherDerived>& other, Scalar prec) const {
    return rotation().isApprox(other.rotation(), prec)
           && ((derived().translationVector()- other.derived().translationVector()).norm() < prec);
};

template<typename Derived>
void TransformBase<Derived>::setIdentity() {
    derived().setIdentity();
}

template<typename Derived>
const Derived TransformBase<Derived>::Identity() {
    return Derived();
}

template<typename Derived>
Derived& TransformBase<Derived>::normalize() {
    derived().normalize();
    return derived();
}

template<typename Derived>
Eigen::Transform<typename TransformBase<Derived>::Scalar, TransformBase<Derived>::Dimension, Eigen::AffineCompact> 
TransformBase<Derived>::transformation() {
    return Eigen::Transform<typename Derived::Scalar, Derived::Dimension, Eigen::AffineCompact>(derived().toRotationMatrix()) *
            Eigen::Transform<typename Derived::Scalar, Derived::Dimension, Eigen::AffineCompact>(derived().translation());
};

/*When you overload a binary operator as a member function of a class the overload is used
 * when the first operand is of the class type.For stream operators, the first operand
 * is the stream and not (usually) the custom class.
*/
template<typename charT, typename traits, typename Derived>
std::basic_ostream<charT,traits>& operator<<(std::basic_ostream<charT,traits>& os, const dcm::details::TransformBase<Derived>& t) {
    os << "Rotation:    " << t.rotation().coeffs()<< std::endl
       << "Translation: " << t.derived().translationVector() <<std::endl;
    return os;
}

}//details
}//DCM



#endif //DCM_TRANSFORMATION
