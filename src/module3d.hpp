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
    typename Kernel::Matrix3	m_rotation;
    typename Kernel::Matrix 	m_diffrot;
    typename Kernel::Quaternion	m_quaternion;
    typename Kernel::Vector3Map	m_normQ;

public:

    void setRotationMap(typename Kernel::Matrix3Map& map, typename Kernel::MatrixMap& diffmap) {
        new(&map) typename Kernel::Matrix3Map(&m_rotation(0,0),3,3,
                                              typename Kernel::DynStride(0,0));
        new(&diffmap) typename Kernel::MatrixMap(&m_diffrot(0,0),3,9,
                typename Kernel::DynStride(0,0));
    };

    typename Kernel::Vector3Map& getNormQuaternionMap() {
        return m_normQ;
    };

    typename Kernel::Quaternion& getQuaternion() {

    };

    void setQuaternionFromNQ() {

    };

    void calculate() {

        typename Kernel::number_type norm = m_normQ.norm(), cn = std::cos(norm), sn = std::sin(norm);

        //typename Kernel::Quaternion Q;
        if(Kernel::isSame(norm, 0)) {

        };
        typename Kernel::Quaternion Q(sn*m_normQ(0)/norm, sn*m_normQ(0)/norm, sn*m_normQ(0)/norm, cn*norm);

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

            MES(Cluster& cl, int par, int eqn) : Sys::Kernel::MappedEquationSystem(par, eqn),
                m_cluster(cl) {};

            virtual void recalculateResidual() {

            }
            virtual void recalculateJacobi() {

            }

            Cluster& m_cluster;
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

            double calculate() {
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
                virtual double calculate(Storage&, Storage&) const = 0;
                virtual placeholder* resetConstraint(Geom first, Geom second) const = 0;
            };

            template< template<typename,typename> class T1, typename T2, typename T3>
            struct holder : public placeholder  {

                holder(const T1<T2,T3> & value)
                    : held(value)   {}

                virtual double calculate(Storage& f, Storage& s) const {
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

        /*
            struct dof_prop {
                typedef edge_property kind;
                typedef Dof<typename Sys::Kernel, Cons> type;
            };*/
        struct vertex_prop {
            typedef Geometry3D kind;
            typedef GlobalVertex type;
        };
        struct edge_prop {
            typedef Constraint3D kind;
            typedef GlobalEdge type;
        };

        typedef mpl::vector<vertex_prop, edge_prop>  properties;

        static void system_init(Sys& sys) {
            sys.m_sheduler.addPreprocessJob(SystemSolver());
        };
    };
};

}

#endif //DCM_GEOMETRY3D_H
