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

#ifndef DCM_GEOMETRY3D_H
#define DCM_GEOMETRY3D_H

#include <boost/mpl/vector.hpp>
#include <boost/mpl/map.hpp>
#include <boost/mpl/less.hpp>
#include <boost/mpl/if.hpp>
#include <boost/mpl/size.hpp>

#include <boost/static_assert.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/preprocessor/iteration/local.hpp>
#include <boost/variant.hpp>

#include <boost/function.hpp>
#include <boost/bind.hpp>

#include "object.hpp"
#include "geometry.hpp"
#include "constraint3d.hpp"
#include "dof.hpp"
#include "clustergraph.hpp"
#include "sheduler.hpp"
#include "traits.hpp"

namespace mpl = boost::mpl;

namespace dcm {

namespace details {

template<typename seq, typename t>
struct distance {
    typedef typename mpl::find<seq, t>::type iterator;
    typedef typename mpl::distance<typename mpl::begin<seq>::type, iterator>::type type;
    BOOST_MPL_ASSERT((mpl::not_< boost::is_same<iterator, typename mpl::end<seq>::type > >));
};

template<typename Sys>
struct ClusterMath {

private:
    typedef typename system_traits<Sys>::Kernel Kernel;
    typedef typename Kernel::number_type Scalar;

    typename Kernel::Matrix3	m_rotation;
    typename Kernel::Matrix39 	m_diffrot;
    typename Kernel::Quaternion	m_quaternion;
    typename Kernel::Vector3Map	m_normQ;

public:
    ClusterMath() : m_normQ(NULL) {};

    void setRotationMap(typename Kernel::Matrix3Map& map, typename Kernel::Matrix39Map& diffmap) {
        new(&map) typename Kernel::Matrix3Map(&m_rotation(0,0),3,3);
        new(&diffmap) typename Kernel::Matrix39Map(&m_diffrot(0,0));
    };

    typename Kernel::Vector3Map& getNormQuaternionMap() {
        return m_normQ;
    };
    typename Kernel::Quaternion& getQuaternion() {
        return m_quaternion;
    };

    void setQuaternionFromNQ() {
        Scalar norm = m_normQ.norm();
        m_quaternion = typename Kernel::Quaternion(m_normQ(0)/norm, m_normQ(1)/norm, m_normQ(2)/norm, norm);
    };

    void calculate() {

        //get the Quaternion for the norm quaternion form and calculate the rotation matrix
        Scalar norm = m_normQ.norm();

        typename Kernel::Quaternion Q(m_normQ(0)/norm, m_normQ(1)/norm, m_normQ(2)/norm, norm);
        if(Kernel::isSame(norm, 0)) {
            Q.setIdentity();
            m_rotation.setIdentity();
            m_diffrot.setZero();
            return;
        };

        Q.normalize();
        m_rotation = Q.toRotationMatrix();

        /* now calculate the gradient quaternions and calculate the diff rotation matrices
         * normQ = (a,b,c)
         * sn0 = ||normQ||^2, n0 = ||normQ||
         *
         * Q = (a/n0, b/n0, c/n0, n0)
         * ||Q|| = n = 1/n0 * sqrt( a^2+b^2+c^2+sn0^2 )
         * n2 = sqrt( a^2+b^2+c^2+sn0^2 )
         *
         * unit Quaternion uQ = (x y z w) = (a/(n0*n), b/(n0*n), c/(n0*n), n0/n)
         * uQ = (a/n2, b/n2, c/n2, sn0/n2)
         */

        const Scalar sn0 = m_normQ.squaredNorm();
        const Scalar n0 = std::sqrt(sn0);
        const Scalar n2 = std::sqrt(sn0 + std::pow(sn0, 2));

        // d(1/n2)/dx and dy and dz
        const Scalar ddn2 = -(1+2*sn0)/std::pow(n2,3);
        const Scalar n2_dda = m_normQ(0)*ddn2;
        const Scalar n2_ddb = m_normQ(1)*ddn2;
        const Scalar n2_ddc = m_normQ(2)*ddn2;

        //dxa = da/dx
        const Scalar dxa = 1/n2 + m_normQ(0)*n2_dda;
        const Scalar dxb = m_normQ(0)*n2_ddb;
        const Scalar dxc = m_normQ(0)*n2_ddc;

        const Scalar dya = m_normQ(1)*n2_dda;
        const Scalar dyb = 1/n2 + m_normQ(1)*n2_ddb;
        const Scalar dyc = m_normQ(1)*n2_ddc;

        const Scalar dza = m_normQ(2)*n2_dda;
        const Scalar dzb = m_normQ(2)*n2_ddb;
        const Scalar dzc = 1/n2 + m_normQ(2)*n2_ddc;

        const Scalar dwa = 2*m_normQ(0)/n2 + sn0*n2_dda;
        const Scalar dwb = 2*m_normQ(1)/n2 + sn0*n2_ddb;
        const Scalar dwc = 2*m_normQ(2)/n2 + sn0*n2_ddc;

        //write in the diffrot matrix, starting with duQ/dx
        m_diffrot(0,0) = -4.0*(Q.y()*dya+Q.z()*dza);
        m_diffrot(0,1) = -2.0*(Q.w()*dza+dwa*Q.z())+2.0*(Q.x()*dya+dxa*Q.y());
        m_diffrot(0,2) = 2.0*(dwa*Q.y()+Q.w()*dya)+2.0*(dxa*Q.z()+Q.x()*dza);
        m_diffrot(1,0) = 2.0*(Q.w()*dza+dwa*Q.z())+2.0*(Q.x()*dya+dxa*Q.y());
        m_diffrot(1,1) = -4.0*(Q.x()*dxa+Q.z()*dza);
        m_diffrot(1,2) = -2.0*(dwa*Q.x()+Q.w()*dxa)+2.0*(dya*Q.z()+Q.y()*dza);
        m_diffrot(2,0) = -2.0*(dwa*Q.y()+Q.w()*dya)+2.0*(dxa*Q.z()+Q.x()*dza);
        m_diffrot(2,1) = 2.0*(dwa*Q.x()+Q.w()*dxa)+2.0*(dya*Q.z()+Q.y()*dza);
        m_diffrot(2,2) = -4.0*(Q.x()*dxa+Q.y()*dya);

        m_diffrot(0,3) = -4.0*(Q.y()*dyb+Q.z()*dzb);
        m_diffrot(0,4) = -2.0*(Q.w()*dzb+dwb*Q.z())+2.0*(Q.x()*dyb+dxb*Q.y());
        m_diffrot(0,5) = 2.0*(dwb*Q.y()+Q.w()*dyb)+2.0*(dxb*Q.z()+Q.x()*dzb);
        m_diffrot(1,3) = 2.0*(Q.w()*dzb+dwb*Q.z())+2.0*(Q.x()*dyb+dxb*Q.y());
        m_diffrot(1,4) = -4.0*(Q.x()*dxb+Q.z()*dzb);
        m_diffrot(1,5) = -2.0*(dwb*Q.x()+Q.w()*dxb)+2.0*(dyb*Q.z()+Q.y()*dzb);
        m_diffrot(2,3) = -2.0*(dwb*Q.y()+Q.w()*dyb)+2.0*(dxb*Q.z()+Q.x()*dzb);
        m_diffrot(2,4) = 2.0*(dwb*Q.x()+Q.w()*dxb)+2.0*(dyb*Q.z()+Q.y()*dzb);
        m_diffrot(2,5) = -4.0*(Q.x()*dxb+Q.y()*dyb);

        m_diffrot(0,6) = -4.0*(Q.y()*dyc+Q.z()*dzc);
        m_diffrot(0,7) = -2.0*(Q.w()*dzc+dwc*Q.z())+2.0*(Q.x()*dyc+dxc*Q.y());
        m_diffrot(0,8) = 2.0*(dwc*Q.y()+Q.w()*dyc)+2.0*(dxc*Q.z()+Q.x()*dzc);
        m_diffrot(1,6) = 2.0*(Q.w()*dzc+dwc*Q.z())+2.0*(Q.x()*dyc+dxc*Q.y());
        m_diffrot(1,7) = -4.0*(Q.x()*dxc+Q.z()*dzc);
        m_diffrot(1,8) = -2.0*(dwc*Q.x()+Q.w()*dxc)+2.0*(dyc*Q.z()+Q.y()*dzc);
        m_diffrot(2,6) = -2.0*(dwc*Q.y()+Q.w()*dyc)+2.0*(dxc*Q.z()+Q.x()*dzc);
        m_diffrot(2,7) = 2.0*(dwc*Q.x()+Q.w()*dxc)+2.0*(dyc*Q.z()+Q.y()*dzc);
        m_diffrot(2,8) = -4.0*(Q.x()*dxc+Q.y()*dyc);

    };
};

}

struct reset {};

template<typename Typelist>
struct Module3D {

    template<typename Sys>
    struct type {
        class Constraint3D;
        class Geometry3D;
        typedef boost::shared_ptr<Geometry3D> Geom;
        typedef mpl::map< mpl::pair<reset, boost::function<void (Geom) > > >  GeomSignal;

        struct MES  : public Sys::Kernel::MappedEquationSystem {

            typedef typename system_traits<Sys>::Cluster Cluster;
	    Cluster& m_cluster;

            MES(Cluster& cl, int par, int eqn) : Sys::Kernel::MappedEquationSystem(par, eqn),
                m_cluster(cl) {};

            virtual void recalculateResidual() {

            }
            virtual void recalculateJacobi() {

            }
        };

        struct SystemSolver : public Job<Sys> {

            typedef typename system_traits<Sys>::Cluster Cluster;

            SystemSolver() {
                Job<Sys>::priority = 1000;
            };

            virtual void execute(Sys& sys, Sheduler<Sys>& shd) {

            };

            void solveCluster(Cluster& cluster) {

                uint parameters, constraints;

                typedef typename boost::graph_traits<Cluster>::vertex_iterator iter;
                std::pair<iter, iter>  it = boost::vertices(cluster);
                for(; it.first != it.second; it.first++) {

                    if(cluster.isCluster(*it.first)) parameters += 3;
                    else parameters += cluster.template getObject<Geometry3D>(*it.first)->m_parameterCount;
                };

                typedef typename boost::graph_traits<Cluster>::edge_iterator e_iter;
                std::pair<e_iter, e_iter>  e_it = boost::edges(cluster);
                for(; e_it.first != e_it.second; e_it.first++)
                    constraints += cluster.getGlobalEdgeCount(*e_it.first);

                MES mes(cluster, parameters, constraints);
            };

        };

        class Geometry3D : public Object<Sys, Geometry3D, GeomSignal > {
            typedef typename boost::make_variant_over< Typelist >::type Variant;
            typedef Object<Sys, Geometry3D, GeomSignal> base;

        public:
            template<typename T>
            Geometry3D(T geometry, Sys& system) : base(system),
                m_geometry(geometry), m_rotation(NULL,0,0), m_parameter(NULL,0), m_diffrot(NULL,0,0)  {};

            template<typename T>
            void set(T geometry) {
                m_geometry = geometry;
                m_parameterCount = geometry_traits<T>::tag::parameters::value;
                base::template emitSignal<reset> (base::shared_from_this());
            };

        protected:
            Variant m_geometry;
            int     m_parameterCount;
            typename Sys::Kernel::Vector    m_original, m_value;
            typename Sys::Kernel::Matrix    m_diffparam;
            typename Sys::Kernel::VectorMap m_parameter;
            typename Sys::Kernel::MatrixMap m_rotation, m_diffrot;

            typename Sys::Kernel::VectorMap& getParameterMap() {
                return m_parameter;
            };
            typename Sys::Kernel::MatrixMap& getRotationMap() {
                return m_rotation;
            };
            typename Sys::Kernel::MatrixMap& getDiffRotationMap() {
                return m_diffrot;
            };

            friend class SystemBuilder;
            friend class Constraint3D;
        };

        //type erasure container for constraints
        class Constraint3D : public Object<Sys, Constraint3D, mpl::map<> > {

            typedef typename system_traits<Sys>::Kernel::number_type Scalar;

        public:
            Constraint3D(Sys& system, Geom f, Geom s) : Object<Sys, Constraint3D, mpl::map<> > (system),
                first(f), second(s), content(0) {

                cf = first->template connectSignal<reset> (boost::bind(&Constraint3D::geometryReset, this, _1));
                cs = second->template connectSignal<reset> (boost::bind(&Constraint3D::geometryReset, this, _1));
            };

            ~Constraint3D()  {
                delete content;
                first->template disconnectSignal<reset>(cf);
                second->template disconnectSignal<reset>(cs);
            }

            template<template<typename, typename> class T>
            void setType() {
                creator<T> creator;
                boost::apply_visitor(creator, first->m_geometry, second->m_geometry);
                content = creator.p;
            };

            Scalar calculate() {
                return content->calculate(first->m_storage, second->m_storage);
            };

        protected:
            void geometryReset(Geom g) {
                placeholder* p = content->resetConstraint(first, second);
                delete content;
                content = p;
            };

            struct placeholder  {

                virtual ~placeholder() {}
                virtual Scalar calculate(Storage&, Storage&) const = 0;
                virtual placeholder* resetConstraint(Geom first, Geom second) const = 0;
            };

            template< template<typename,typename> class T1, typename T2, typename T3>
            struct holder : public placeholder  {

                holder(const T1<T2,T3> & value)
                    : held(value)   {}

                virtual Scalar calculate(Storage& f, Storage& s) const {
                    return held.calculate(f,s);
                };

                virtual placeholder* resetConstraint(Geom first, Geom second) const {
                    creator<T1> creator;
                    boost::apply_visitor(creator, first->m_geometry, second->m_geometry);
                    return creator.p;
                };

                T1<T2,T3>  held;
            };

            template<template<typename, typename> class T>
            struct creator : public boost::static_visitor<void> {

                template<typename T1, typename T2>
                void operator()(const T1&, const T2&) {
                    typedef T<typename geometry_traits<T1>::tag, typename geometry_traits<T2>::tag> type;
                    p = new holder< T, typename geometry_traits<T1>::tag, typename geometry_traits<T2>::tag > (type());
                };
                placeholder* p;
            };

            placeholder* content;
            Geom first, second;
            Connection cf, cs;

        };

        typedef boost::shared_ptr<Constraint3D> Cons;
        typedef mpl::vector<Geometry3D, Constraint3D> objects;

        struct inheriter {
            inheriter() {
                m_this = ((Sys*) this);
            };

            template<typename T>
            Geom createGeometry3D(T geom) {

                Geom g(new Geometry3D(geom, * ((Sys*) this)));
                fusion::vector<LocalVertex, GlobalVertex> res = m_this->m_cluster.addVertex();
                m_this->m_cluster.template setObject<Geometry3D> (fusion::at_c<0> (res), g);
                g->template setProperty<vertex_prop>(fusion::at_c<1>(res));
                return g;
            };

            template<template<typename, typename> class T>
            Cons createConstraint3D(Geom first, Geom second) {

                Cons c(new Constraint3D(* ((Sys*) this), first, second));
                c->template setType<T>();
                fusion::vector<LocalEdge, GlobalEdge, bool, bool> res;
                res = m_this->m_cluster.addEdge(first->template getProperty<vertex_prop>(),
                                                second->template getProperty<vertex_prop>());
                c->template setProperty<edge_prop>(fusion::at_c<1>(res));
                return c;
            };

        private:
            Sys* m_this;
        };


        struct math_prop {
            typedef cluster_property kind;
            typedef details::ClusterMath<Sys> type;
        };
        struct vertex_prop {
            typedef Geometry3D kind;
            typedef GlobalVertex type;
        };
        struct edge_prop {
            typedef Constraint3D kind;
            typedef GlobalEdge type;
        };

        typedef mpl::vector<vertex_prop, edge_prop, math_prop>  properties;

        static void system_init(Sys& sys) {
            sys.m_sheduler.addProcessJob(SystemSolver());
        };
    };
};

}

#endif //DCM_GEOMETRY3D_H
