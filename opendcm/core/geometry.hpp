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


#ifndef GCM_GEOMETRY_H
#define GCM_GEOMETRY_H

#include <iostream>

#include <eigen3/Eigen/Core>

#include <boost/type_traits.hpp>
#include <boost/mpl/assert.hpp>
#include <boost/mpl/fold.hpp>
#include <boost/mpl/set.hpp>
#include <boost/mpl/less.hpp>
#include <boost/mpl/or.hpp>
#include <boost/mpl/not.hpp>
#include <boost/mpl/bool.hpp>

#include <boost/fusion/include/as_vector.hpp>
#include <boost/fusion/include/mpl.hpp>
#include <boost/concept_check.hpp>
#include <boost/graph/graph_concepts.hpp>

#include <boost/variant.hpp>

#include "object.hpp"
#include "traits.hpp"
#include "logging.hpp"


namespace mpl = boost::mpl;
namespace fusion = boost::fusion;

namespace dcm {

namespace tag {

struct undefined {
    typedef mpl::int_<0> parameters;
    typedef mpl::int_<0> transformations;
};

//we need to order tags, this values make it easy for module tags
namespace weight {
struct point : mpl::int_<1> {};
struct line  : mpl::int_<2> {};
struct plane : mpl::int_<3> {};
struct cylinder : mpl::int_<4> {};
}
}

struct orderd_bracket_accessor {

    template<typename Scalar, int ID, typename T>
    Scalar get(T& t) {
        return t[ID];
    };
    template<typename Scalar, int ID, typename T>
    void set(Scalar value, T& t) {
        t[ID] = value;
    };
};

struct orderd_roundbracket_accessor {

    template<typename Scalar, int ID, typename T>
    Scalar get(T& t) {
        return t(ID);
    };
    template<typename Scalar, int ID,  typename T>
    void set(Scalar value, T& t) {
        t(ID) = value;
    };
};

//tag ordering
template<typename T1, typename T2>
struct tag_order {

    BOOST_MPL_ASSERT((mpl::not_< mpl::or_<
                      boost::is_same< typename T1::weight, mpl::int_<0> >,
                      boost::is_same< typename T2::weight, mpl::int_<0> >  >  >));

    typedef typename mpl::less<typename T2::weight, typename T1::weight>::type swapt;
    typedef typename mpl::if_<swapt, T2, T1>::type first_tag;
    typedef typename mpl::if_<swapt, T1, T2>::type second_tag;
};


//template<typename T1, typename T2>
//struct type_order : public tag_order< typename geometry_traits<T1>::tag, typename geometry_traits<T2>::tag > {};


template< typename T>
struct geometry_traits {
    BOOST_MPL_ASSERT_MSG(false, NO_GEOMETRY_TRAITS_SPECIFIED_FOR_TYPE, (T));
};

struct reset {}; 	//signal namespace

namespace detail {

template< typename Sys, typename Derived, typename GeometrieTypeList, typename Signals, int Dimension>
class Geometry : public Object<Sys, Derived, Signals > {

    typedef typename boost::make_variant_over< GeometrieTypeList >::type Variant;
    typedef Object<Sys, Derived, Signals> 		Base;
    typedef typename system_traits<Sys>::Kernel 	Kernel;
    typedef typename system_traits<Sys>::Cluster 	Cluster;
    typedef typename Kernel::number_type 		Scalar;
    typedef typename Kernel::DynStride 			DS;
    typedef typename Kernel::Transform3D		Transform;
    typedef typename Kernel::DiffTransform3D		DiffTransform;

#ifdef USE_LOGGING
protected:
    src::logger log;
#endif

public:
    template<typename T>
    Geometry(T geometry, Sys& system) : Base(system), m_isInCluster(false),
        m_geometry(geometry),  m_parameter(NULL,0,DS(0,0)),
        m_clusterFixed(false), m_init(false) {

#ifdef USE_LOGGING
        log.add_attribute("Tag", attrs::constant< std::string >("Geometry3D"));
#endif

        init<T>(geometry);
    };

    template<typename T>
    void set(T geometry) {
        m_geometry = geometry;
        init<T>(geometry);
        Base::template emitSignal<reset> (Base::shared_from_this());
    };

    template<typename Visitor>
    typename Visitor::result_type apply(Visitor& vis) {
        return boost::apply_visitor(vis, m_geometry);
    };

public:
    Variant m_geometry; //Variant holding the real geometry type
    int     m_BaseParameterCount; //count of the parameters the variant geometry type needs
    int     m_parameterCount; //count of the used parameters (when in cluster:6, else m_BaseParameterCount)
    int     m_offset; //the starting point of our parameters in the math system parameter vector
    int     m_rotations; //count of rotations to be done when original vector gets rotated
    int     m_translations; //count of translations to be done when original vector gets rotated
    bool    m_isInCluster, m_clusterFixed, m_init;
    typename Sys::Kernel::Vector      m_toplocal; //the local value in the toplevel cluster used for cuttent solving
    typename Sys::Kernel::Vector      m_global; //the global value outside of all clusters
    typename Sys::Kernel::Vector      m_rotated; //the global value as the rotation of toplocal (used as temp)
    typename Sys::Kernel::Matrix      m_diffparam; //gradient vectors combined as matrix when in cluster
    typename Sys::Kernel::VectorMap   m_parameter; //map to the parameters in the solver

    template<typename T>
    void init(T& t) {
        m_BaseParameterCount = geometry_traits<T>::tag::parameters::value;
        m_parameterCount = m_BaseParameterCount;
        m_rotations = geometry_traits<T>::tag::rotations::value;
        m_translations = geometry_traits<T>::tag::translations::value;

        m_toplocal.setZero(m_parameterCount);
        m_global.resize(m_parameterCount);
        m_rotated.resize(m_parameterCount);
	m_rotated.setZero();

        m_diffparam.resize(m_parameterCount,6);
        m_diffparam.setZero();

        (typename geometry_traits<T>::modell()).template extract<Scalar,
        typename geometry_traits<T>::accessor >(t, m_global);

        //new value which is not set into parameter, so init is false
        m_init = false;

#ifdef USE_LOGGING
        BOOST_LOG(log) << "Init: "<<m_global.transpose();
#endif

    }

    typename Sys::Kernel::VectorMap& getParameterMap() {
        m_isInCluster = false;
        m_parameterCount = m_BaseParameterCount;
        return m_parameter;
    }
    void initMap() {
        //when direct parameter solving the global value is wanted (as it's the initial rotation*toplocal)
        m_parameter = m_global;
        m_init = true;
    };

    void setClusterMode(bool iscluster, bool isFixed) {
        m_isInCluster = iscluster;
        m_clusterFixed = isFixed;
        if(iscluster) {
            //we are in cluster, therfore the parameter map should not point to a solver value but to
            //the rotated original value;
            new(&m_parameter) typename Sys::Kernel::VectorMap(&m_rotated(0), m_parameterCount, DS(1,1));
            //the local value is the global one as no transformation was applied  yet
            m_toplocal = m_global;
        } else new(&m_parameter) typename Sys::Kernel::VectorMap(&m_global(0), m_parameterCount, DS(1,1));

    }
    bool getClusterMode() {
        return m_isInCluster;
    };
    bool isClusterFixed() {
        return m_clusterFixed;
    };

    void recalculate(DiffTransform& trans) {
        if(!m_isInCluster) return;

        for(int i=0; i!=m_rotations; i++) {
            //first rotate the original to the transformed value
            m_rotated.block(i*Dimension,0,Dimension,1) = trans.rotation()*m_toplocal.template segment<Dimension>(i*Dimension);

            //now calculate the gradient vectors and add them to diffparam
            for(int j=0; j<Dimension; j++)
                m_diffparam.block(i*Dimension,j,Dimension,1) = trans.differential().block(0,j*3,Dimension,Dimension) * m_toplocal.template segment<Dimension>(i*Dimension);
        }
        //after rotating the needed parameters we translate the stuff that needs to be moved
        for(int i=0; i!=m_translations; i++) {
            m_rotated.block(i*Dimension,0,Dimension,1) += trans.translation().vector();
            m_rotated.block(i*Dimension,0,Dimension,1) *= trans.scaling();
            //calculate the gradient vectors and add them to diffparam
            m_diffparam.block(i*Dimension,Dimension,Dimension,Dimension).setIdentity();
        }
        //basicly no needed as this case appears only from within finishCalculation
        //and therefore the diffparam is not needed anymore
        //if(trans.scaling() != 1.) m_diffparam *= trans.scaling();
    }

    typename Kernel::Vector3 getPoint() {
        return m_toplocal.template segment<Dimension>(0);
    }

    //visitor to write the calculated value into the variant
    struct apply_visitor : public boost::static_visitor<void> {

        apply_visitor(typename Kernel::Vector& v) : value(v) {};
        template <typename T>
        void operator()(T& t) const  {
            (typename geometry_traits<T>::modell()).template inject<typename Kernel::number_type,
            typename geometry_traits<T>::accessor >(t, value);
        }
        typename Kernel::Vector& value;
    };

    //use m_value or parametermap as new value, dependend on the solving mode
    void finishCalculation() {
        //if fixed nothing needs to be changed
        if(m_isInCluster) {
            //recalculate(1.); //remove scaling to get right global value
            m_global = m_rotated;
        }
        //TODO:non cluster paramter scaling
        else m_global = m_parameter;
        apply_visitor v(m_global);
        apply(v);
    };


    //normal transformation
    void transform(const Transform& t) {

        if(m_isInCluster)
            transform(t, m_toplocal);
        else if(m_init)
            transform(t, m_rotated);
        else
            transform(t, m_global);
    }

    template<typename VectorType>
    void transform(const Transform& t, VectorType& vec) {

        //everything that needs to be translated needs to be fully transformed
        for(int i=0; i!=m_translations; i++) {
            typename Kernel::Vector3 v = vec.template segment<Dimension>(i*Dimension);
            vec.template segment<Dimension>(i*Dimension) = t*v;
        }

        for(int i=m_translations; i!=m_rotations; i++) {
            typename Kernel::Vector3 v = vec.template segment<Dimension>(i*Dimension);
            vec.template segment<Dimension>(i*Dimension) = t.rotate(v);
        }

#ifdef USE_LOGGING
        BOOST_LOG(log) << "Transformed with cluster: "<<m_isInCluster
                       << ", init: "<<m_init<<" into: "<< vec.transpose();
#endif
    }

    void scale(Scalar value) {

        for(int i=0; i!=m_translations; i++)
            m_parameter.template segment<Dimension>(i*Dimension) *= 1./value;

    };

};

}
}

#endif // GCM_GEOMETRY_H
