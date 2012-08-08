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
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/

#ifndef GCM_MODULE_PART_H
#define GCM_MODULE_PART_H

#include "opengcm/Core"
#include "opengcm/core/traits.hpp"
#include "opengcm/core/clustergraph.hpp"
#include "opengcm/core/property.hpp"
#include "opengcm/Module3D"

#include <boost/mpl/assert.hpp>
#include <boost/utility/enable_if.hpp>

namespace mpl = boost::mpl;

namespace gcm {

enum { clusterPart = 110};

enum CoordinateFrame {Local, Global};

template<typename Typelist, typename Identifier = No_Identifier>
struct ModulePart {

    template<typename Sys>
    struct type {


        class Part;
        struct PrepareCluster;
        struct EvaljuateCluster;
        typedef boost::shared_ptr<Part> Partptr;
        typedef mpl::map< >  PartSignal;

        class Part_base : public Object<Sys, Part, PartSignal > {
        protected:
            using Object<Sys, Part, PartSignal >::m_system;

            //check if we have module3d in this system
            typedef typename system_traits<Sys>::template getModule<m3d> getM;
            typedef typename system_traits<Sys>::template getModule<m3d>::type module3d;
            BOOST_MPL_ASSERT((mpl::not_<boost::is_same<module3d, mpl::void_> >));

            //define what we need
            typedef typename module3d::Geometry3D Geometry3D;
            typedef boost::shared_ptr<Geometry3D> Geom;

            typedef typename boost::make_variant_over< Typelist >::type Variant;
            typedef Object<Sys, Part, PartSignal> base;
            typedef typename system_traits<Sys>::Kernel Kernel;
            typedef typename system_traits<Sys>::Cluster Cluster;
            typedef typename Kernel::number_type Scalar;

        public:
            template<typename T>
            Part_base(T geometry, Sys& system, Cluster& cluster) : base(system),
                m_geometry(geometry), m_cluster(cluster)  {

                (typename geometry_traits<T>::modell()).template extract<typename Part_base::Scalar,
                typename geometry_traits<T>::accessor >(geometry, Part_base::m_quaternion, Part_base::m_translation);

                Part_base::m_quaternion.normalize();		
		//the cluster needs initial values but they are set by preprocess job
            };

            template<typename Visitor>
            typename Visitor::result_type apply(Visitor& vis) {
                return boost::apply_visitor(vis, m_geometry);
            };

        public:
            Variant m_geometry;
            typename Kernel::Quaternion m_quaternion;
            typename Kernel::Vector3 m_translation;
            Cluster& m_cluster;

            template<typename T>
            Geom addGeometry(T geom, CoordinateFrame frame = Global) {
                Geom g(new Geometry3D(geom, m_system));
                if(frame == Local)
                    g->transformGlobal(m_quaternion.toRotationMatrix(), m_translation);

                fusion::vector<LocalVertex, GlobalVertex> res = m_cluster.addVertex();
                m_cluster.template setObject<Geometry3D> (fusion::at_c<0> (res), g);
                g->template setProperty<typename module3d::vertex_prop>(fusion::at_c<1>(res));
                m_system.template objectVector<Geometry3D>().push_back(g);

                return g;
            };

            //visitor to write the calculated value into the variant
            struct apply_visitor : public boost::static_visitor<void> {

                apply_visitor(typename Kernel::Quaternion& v1,typename Kernel::Vector3& v2) : quaternion(v1),
                    translation(v2) {};
                template <typename T>
                void operator()(T& t) const  {
                    (typename geometry_traits<T>::modell()).template inject<typename Part_base::Scalar,
                    typename geometry_traits<T>::accessor >(t, quaternion, translation);
                }
                typename Kernel::Vector3& translation;
                typename Kernel::Quaternion& quaternion;
            };

            void finishCalculation() {
                m_quaternion.normalize();
                apply_visitor vis(m_quaternion, m_translation);
                apply(vis);
            };
        };

        class Part_noid : public Part_base {

        public:
            template<typename T>
            Part_noid(T geometry, Sys& system, typename Part_base::Cluster& cluster) : Part_base(geometry, system, cluster) {};

            template<typename T>
            typename Part_base::Geom addGeometry3D(T geom, CoordinateFrame frame = Global) {
                return  Part_base::addGeometry(geom, frame);
            };

            template<typename T>
            void set(T geometry) {
                Part_base::m_geometry = geometry;
                (typename geometry_traits<T>::modell()).template extract<typename Part_base::Scalar,
                typename geometry_traits<T>::accessor >(geometry, Part_base::m_quaternion, Part_base::m_translation);
            };
        };
        class Part_id : public Part_base {

            Identifier m_id;
        public:
            template<typename T>
            Part_id(T geometry, Sys& system,  typename Part_base::Cluster& cluster) : Part_base(geometry, system, cluster) {};

            template<typename T>
            typename Part_base::Geom addGeometry3D(T geom, Identifier id, CoordinateFrame frame = Local) {

                typename Part_base::Geom g = Part_base::addGeometry(geom, frame);
                g->setIdentifier(id);
                return g;
            };

            template<typename T>
            void set(T geometry, Identifier id) {
                Part_base::m_geometry = geometry;
                m_id = id;
                (typename geometry_traits<T>::modell()).template extract< typename Part_base::Scalar,
                typename geometry_traits<T>::accessor >(geometry, Part_base::m_quaternion, Part_base::m_translation);
            };

            bool hasGeometry3D(Identifier id) {
                typename Part_base::Geom g = Part_base::m_system.getGeometry3D(id);
                if(!g) return false;

                //get the global vertex and check if it is a child of the part cluster
                GlobalVertex v = g->template getProperty<typename Part_base::module3d::vertex_prop>();
                return Part_base::m_cluster.getLocalVertex(v).second;
            };

            Identifier& getIdentifier() {
                return m_id;
            };

            void setIdentifier(Identifier id) {
                m_id = id;
            };
        };

        struct Part : public mpl::if_<boost::is_same<Identifier, No_Identifier>, Part_noid, Part_id>::type {

            typedef typename mpl::if_<boost::is_same<Identifier, No_Identifier>, Part_noid, Part_id>::type base;
            template<typename T>
            Part(T geometry, Sys& system, typename base::Cluster& cluster) : base(geometry, system, cluster) {};

            friend struct PrepareCluster;
            friend struct EvaljuateCluster;
        };


        struct inheriter_base {

            inheriter_base() {
                m_this = ((Sys*) this);
            };

        protected:
            Sys* m_this;

            template<typename T>
            Partptr partCreation(T geometry) {
                typedef typename system_traits<Sys>::Cluster Cluster;
                std::pair<Cluster&, LocalVertex>  res = m_this->m_cluster.createCluster();
                Partptr p(new Part(geometry, * ((Sys*) this), res.first));

                m_this->m_cluster.template setObject<Part> (res.second, p);
                m_this->template objectVector<Part>().push_back(p);

                res.first.template setClusterProperty<type_prop>(clusterPart);
                return p;
            }
        };

        struct inheriter_noid : public inheriter_base {
            template<typename T>
            Partptr createPart(T geometry) {
                return inheriter_base::partCreation(geometry);
            };
        };

        struct inheriter_id : public inheriter_base {

            template<typename T>
            Partptr createPart(T geometry, Identifier id) {
                Partptr p = inheriter_base::partCreation(geometry);
                p->setIdentifier(id);
                return p;
            };

            bool hasPart(Identifier id) {
                if(getPart(id)) return true;
                return false;
            };

            Partptr getPart(Identifier id) {
                std::vector< Partptr >& vec = inheriter_base::m_this->template objectVector<Part>();
                typedef typename std::vector<Partptr>::iterator iter;
                for(iter it=vec.begin(); it!=vec.end(); it++) {
                    if(compare_traits<Identifier>::compare((*it)->getIdentifier(), id)) return *it;
                };
                return Partptr();
            };
        };

        struct inheriter : public mpl::if_<boost::is_same<Identifier, No_Identifier>, inheriter_noid, inheriter_id>::type {};

        typedef mpl::vector<>  properties;
        typedef mpl::vector<Part>  objects;

        struct PrepareCluster : public Job<Sys> {

            typedef typename system_traits<Sys>::Cluster Cluster;
            typedef typename system_traits<Sys>::Kernel Kernel;
            typedef typename system_traits<Sys>::template getModule<m3d>::type module3d;

            PrepareCluster() {
                Job<Sys>::priority = 1000;
            };

            virtual void execute(Sys& sys) {
                //get all parts and set their values to the cluster's
                typedef typename std::vector<Partptr>::iterator iter;
                for(iter it = sys.template begin<Part>(); it != sys.template end<Part>(); it++) {

                    details::ClusterMath<Sys>& cm = (*it)->m_cluster.template getClusterProperty<typename module3d::math_prop>();
                    cm.getQuaternion() = (*it)->m_quaternion;
                    cm.getTranslation() = (*it)->m_translation;
                };
            };
        };

        struct EvaljuateCluster : public Job<Sys> {

            typedef typename system_traits<Sys>::Cluster Cluster;
            typedef typename system_traits<Sys>::Kernel Kernel;
            typedef typename system_traits<Sys>::template getModule<m3d>::type module3d;

            EvaljuateCluster() {
                Job<Sys>::priority = 1000;
            };

            virtual void execute(Sys& sys) {
                //get all parts and set their values to the cluster's
                typedef typename std::vector<Partptr>::iterator iter;
                for(iter it = sys.template begin<Part>(); it != sys.template end<Part>(); it++) {

                    details::ClusterMath<Sys>& cm = (*it)->m_cluster.template getClusterProperty<typename module3d::math_prop>();
                    (*it)->m_quaternion = cm.getQuaternion();
                    (*it)->m_translation = cm.getTranslation();
                    (*it)->finishCalculation();		    
                };
            };
        };

        static void system_init(Sys& sys) {
            sys.m_sheduler.addPreprocessJob(new PrepareCluster());
            sys.m_sheduler.addPostprocessJob(new EvaljuateCluster());
        };
    };
};

}

#endif //GCM_MODULEPART_H




