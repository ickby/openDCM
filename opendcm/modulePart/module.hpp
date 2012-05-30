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

#ifndef DCM_MODULE_PART_H
#define DCM_MODULE_PART_H

#include "opendcm/Core"
#include "opendcm/core/traits.hpp"
#include "opendcm/core/clustergraph.hpp"
#include "opendcm/core/property.hpp"
#include "opendcm/Module3D"

#include <boost/mpl/assert.hpp>

namespace mpl = boost::mpl;

namespace dcm {

enum { clusterPart = 110};

template<typename Typelist>
struct ModulePart {

    template<typename Sys>
    struct type {




        class Part;
        typedef boost::shared_ptr<Part> Partptr;
        typedef mpl::map< >  PartSignal;

        class Part : public Object<Sys, Part, PartSignal > {

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
            Part(T geometry, Sys& system, Cluster& cluster) : base(system),
                m_geometry(geometry), m_cluster(cluster)  {};

            template<typename T>
            Geom addGeometry3D(T geom) {

                Geom g(new Geometry3D(geom, * ((Sys*) this)));
                fusion::vector<LocalVertex, GlobalVertex> res = m_cluster.addVertex();
                m_cluster.template setObject<Geometry3D> (fusion::at_c<0> (res), g);
                g->template setProperty<typename module3d::vertex_prop>(fusion::at_c<1>(res));
                //((Sys*) this)->template objectVector<Geometry3D>().push_back(g);
                return g;
            };

            template<typename T>
            void set(T geometry) {
                m_geometry = geometry;
                (typename geometry_traits<T>::modell()).template extract<Scalar,
                typename geometry_traits<T>::accessor >(geometry, m_quaternion);
            };

            template<typename Visitor>
            typename Visitor::result_type apply(Visitor& vis) {
                return boost::apply_visitor(vis, m_geometry);
            };



        private:
            Variant m_geometry;
            typename Kernel::Quaternion m_quaternion;
            Cluster& m_cluster;
        };

        struct inheriter {

            inheriter() {
                m_this = ((Sys*) this);
            };


            template<typename T>
            Partptr createPart(T geometry) {

                typedef typename system_traits<Sys>::Cluster Cluster;
                std::pair<Cluster&, LocalVertex>  res = m_this->m_cluster.createCluster();
                Partptr p(new Part(geometry, * ((Sys*) this), res.first));

                m_this->m_cluster.template setObject<Part> (res.second, p);
                m_this->template objectVector<Part>().push_back(p);

                res.first.template setClusterProperty<type_prop>(clusterPart);
                return p;
            };

        private:
            Sys* m_this;
        };

        typedef mpl::vector<>  properties;
        typedef mpl::vector<Part>  objects;

        static void system_init(Sys& sys) {
            //sys.m_sheduler.addProcessJob(new SystemSolver());
        };
    };
};

}

#endif //DCM_MODULEPART_H



