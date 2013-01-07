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

#ifndef GCM_MODULE_3D_H
#define GCM_MODULE_3D_H

#include <boost/mpl/vector.hpp>
#include <boost/mpl/map.hpp>
#include <boost/mpl/less.hpp>
#include <boost/mpl/if.hpp>
#include <boost/mpl/size.hpp>

#include <boost/fusion/include/as_vector.hpp>
#include <boost/fusion/include/mpl.hpp>
#include <boost/fusion/include/at.hpp>

#include <boost/static_assert.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/preprocessor/iteration/local.hpp>
#include <boost/variant.hpp>

#include <boost/function.hpp>
#include <boost/bind.hpp>

#include "opendcm/Core"
#include "opendcm/core/object.hpp"
#include "opendcm/core/clustergraph.hpp"
#include "opendcm/core/sheduler.hpp"
#include "opendcm/core/traits.hpp"
#include "opendcm/core/geometry.hpp"
#include "geometry.hpp"
#include "distance.hpp"
#include "parallel.hpp"
#include "angle.hpp"
#include "dof.hpp"

static int counter = 0;

namespace mpl = boost::mpl;

namespace dcm {

namespace details {

enum { cluster3D = 100};

template<typename seq, typename t>
struct distance {
    typedef typename mpl::find<seq, t>::type iterator;
    typedef typename mpl::distance<typename mpl::begin<seq>::type, iterator>::type type;
    BOOST_MPL_ASSERT((mpl::not_< boost::is_same<iterator, typename mpl::end<seq>::type > >));
};
}

struct m3d {}; 	//base of module3d::type to allow other modules check for it


}//dcm

//needs to be here to access m3d struct
#include "clustermath.hpp"

namespace dcm {

template<typename Typelist, typename Identifier = No_Identifier>
struct Module3D {

    template<typename Sys>
    struct type : m3d {
        class Constraint3D;
        class Geometry3D;
        typedef boost::shared_ptr<Geometry3D> Geom;
        typedef boost::shared_ptr<Constraint3D> Cons;

        typedef mpl::map< mpl::pair<reset, boost::function<void (Geom) > > >  GeomSignal;
        typedef mpl::map<  >  ConsSignal;

        struct MES  : public system_traits<Sys>::Kernel::MappedEquationSystem {

	    typedef typename system_traits<Sys>::Kernel::MappedEquationSystem Base;
            typedef typename system_traits<Sys>::Cluster Cluster;
            Cluster& m_cluster;

            MES(Cluster& cl, int par, int eqn) : Base(par, eqn),
                  m_cluster(cl) {};

            virtual void recalculate() {

                //first calculate all clusters
                typedef typename Cluster::cluster_iterator citer;
                std::pair<citer, citer> cit = m_cluster.clusters();
                for(; cit.first != cit.second; cit.first++) {

                    if(!(*cit.first).second->template getClusterProperty<fix_prop>()) 
                        (*cit.first).second->template getClusterProperty<math_prop>().recalculate();

                };

                //with everything updated just nicely we can compute the constraints
                typedef typename Cluster::template object_iterator<Constraint3D> oiter;
                typedef typename boost::graph_traits<Cluster>::edge_iterator eiter;
                std::pair<eiter, eiter>  eit = boost::edges(m_cluster);
                for(; eit.first != eit.second; eit.first++) {
                    //as always: every local edge can hold multiple global ones, so iterate over all constraints
                    //hold by the individual edge
                    std::pair< oiter, oiter > oit = m_cluster.template getObjects<Constraint3D>(*eit.first);
                    for(; oit.first != oit.second; oit.first++) {
                        if(*oit.first)
                            (*oit.first)->calculate(Base::Scaling);
                    }
                }
            }
        };

        struct SystemSolver : public Job<Sys> {

            typedef typename system_traits<Sys>::Cluster Cluster;
            typedef typename system_traits<Sys>::Kernel Kernel;
            typedef typename Kernel::number_type Scalar;

            SystemSolver() {
                Job<Sys>::priority = 1000;
            };

            virtual void execute(Sys& sys) {
                solveCluster(sys.m_cluster, sys);
            };

            void solveCluster(Cluster& cluster, Sys& sys) {

                //set out and solve all relevant subclusters
                typedef typename Cluster::cluster_iterator citer;
                std::pair<citer, citer> cit = cluster.clusters();
                for(; cit.first != cit.second; cit.first++) {

                    Cluster& c = *(*cit.first).second;
                    if(c.template getClusterProperty<changed_prop>() &&
                            c.template getClusterProperty<type_prop>() == details::cluster3D)
                        solveCluster(c, sys);
                }

                int params=0, constraints=0;
                typename Kernel::number_type scale = 1;

                //get the ammount of parameters and constraint equations we need
                typedef typename boost::graph_traits<Cluster>::vertex_iterator iter;
                std::pair<iter, iter>  it = boost::vertices(cluster);
                for(; it.first != it.second; it.first++) {

                    //when cluster and not fixed it has trans and rot parameter
                    if(cluster.isCluster(*it.first)) {
                        if(!cluster.template getSubclusterProperty<fix_prop>(*it.first)) {
                            params += 6;
                        }
                    } else {
                        params += cluster.template getObject<Geometry3D>(*it.first)->m_parameterCount;
                    };
                }

                //count the equations in the constraints
                typedef typename Cluster::template object_iterator<Constraint3D> ocit;
                typedef typename boost::graph_traits<Cluster>::edge_iterator e_iter;
                std::pair<e_iter, e_iter>  e_it = boost::edges(cluster);
                for(; e_it.first != e_it.second; e_it.first++) {
                    std::pair< ocit, ocit > it = cluster.template getObjects<Constraint3D>(*e_it.first);
                    for(; it.first != it.second; it.first++)
                        constraints += (*it.first)->equationCount();
                };


                //initialise the system with now known size
                //std::cout<<"constraints: "<<constraints<<", params: "<<params+rot_params+trans_params<<std::endl;
                MES mes(cluster, params, constraints);

                //iterate all geometrys again and set the needed maps
                it = boost::vertices(cluster);
                for(; it.first != it.second; it.first++) {

                    if(cluster.isCluster(*it.first)) {
                        Cluster& c = cluster.getVertexCluster(*it.first);
                        details::ClusterMath<Sys>& cm =  c.template getClusterProperty<math_prop>();
                        //only get maps and propagate downstream if not fixed
                        if(!c.template getClusterProperty<fix_prop>()) {
                            //set norm Quaternion as map to the parameter vector
                            int offset = mes.setParameterMap(cm.getNormQuaternionMap());
                            //set translation as map to the parameter vector
                            mes.setParameterMap(cm.getTranslationMap());
                            //write initail values to the parameter maps
                            //remember the parameter offset as all downstream geometry must use this offset
                            cm.setParameterOffset(offset);
                            //wirte initial values
                            cm.initMaps();
                        } else cm.initFixMaps();

                        //map all geometrie within that cluster to it's rotation matrix
                        //for collecting all geometries which need updates
                        cm.clearGeometry();

                        //to allow a corect calculation of geometries toplocal value we need the quaternion
                        //which transforms from toplevel to this "to be solved" cluster and the aquivalent translation
                       	typename Kernel::Transform3D trans;
                        cm.mapClusterDownstreamGeometry(c, trans, cm);


                    } else {
                        Geom g = cluster.template getObject<Geometry3D>(*it.first);
                        int offset = mes.setParameterMap(g->m_parameterCount, g->getParameterMap());
                        g->m_offset = offset;
                        //init the parametermap with initial values
                        g->initMap();
                    }
                }

                //and now the constraints to set the residual and gradient maps
                typedef typename Cluster::template object_iterator<Constraint3D> oiter;
                e_it = boost::edges(cluster);
                for(; e_it.first != e_it.second; e_it.first++) {


                    //as always: every local edge can hold multiple global ones, so iterate over all constraints
                    //hold by the individual edge
                    std::pair< oiter, oiter > oit = cluster.template getObjects<Constraint3D>(*e_it.first);
                    for(; oit.first != oit.second; oit.first++) {

                        //set the maps
                        Cons c = *oit.first;
                        if(c) c->setMaps(mes);
                        //TODO: else throw (as every global edge was counted as one equation)
                    }
                }

                //get the maximal scale
                Scalar sc = 0;
                for(cit = cluster.clusters(); cit.first != cit.second; cit.first++) {
                    //fixed cluster are irrelevant for scaling
                    if((*cit.first).second->template getClusterProperty<fix_prop>()) continue;

                    //get the biggest scale factor
                    const Scalar s = (*cit.first).second->template getClusterProperty<math_prop>().calculateClusterScale();
                    sc = (s>sc) ? s : sc;
                }
                //if no scaling-value returned we can use 1
                sc = (Kernel::isSame(sc,0)) ? 1. : sc;
                //std::cout<<"scaling: "<<sc<<std::endl;
                //Base::Console().Message("Scale is %f\n", sc);
                it = boost::vertices(cluster);
                for(; it.first != it.second; it.first++) {

                    if(cluster.isCluster(*it.first)) {
                        Cluster& c = cluster.getVertexCluster(*it.first);
                        c.template getClusterProperty<math_prop>().applyClusterScale(sc,
                                c.template getClusterProperty<fix_prop>());
                    }
                    else {
                        Geom g = cluster.template getObject<Geometry3D>(*it.first);
                        g->scale(sc);
                    }
                }
                mes.Scaling = 1./sc;
                //std::cout<<"if"<<std::endl;


                //now it's time to solve
                Kernel::solve(mes);

                //std::cout<<"Residual after solving: "<<mes.Residual.norm()<<std::endl;
                //std::cout<<"mes scaling: "<<mes.Scaling<<std::endl;

                //now go to all relevant geometries and clusters and write the values back
                it = boost::vertices(cluster);
                for(; it.first != it.second; it.first++) {

                    if(cluster.isCluster(*it.first)) {
                        Cluster& c = cluster.getVertexCluster(*it.first);
                        if(!cluster.template getSubclusterProperty<fix_prop>(*it.first))
                            c.template getClusterProperty<math_prop>().finishCalculation();
			else
			    c.template getClusterProperty<math_prop>().finishFixCalculation();

                        std::vector<Geom>& vec = c.template getClusterProperty<math_prop>().getGeometry();
                        for(typename std::vector<Geom>::iterator vit = vec.begin(); vit != vec.end(); vit++)
                            (*vit)->finishCalculation();

                    } else {
                        Geom g = cluster.template getObject<Geometry3D>(*it.first);
			g->scale(1./sc);
                        g->finishCalculation();
                    }
                }

                //we have solved this cluster
                cluster.template setClusterProperty<changed_prop>(false);

            };

        };

        template<typename Derived>
        class Geometry3D_id : public detail::Geometry<Sys, Derived, Typelist, GeomSignal, 3> {

            typedef detail::Geometry<Sys, Derived, Typelist, GeomSignal, 3> Base;

            Identifier m_id;
#ifdef USE_LOGGING
            attrs::mutable_constant< std::string > log_id;
#endif
        public:
            template<typename T>
            Geometry3D_id(T geometry, Sys& system) : Base(geometry, system)
#ifdef USE_LOGGING
                , log_id("No ID")
#endif
            {

#ifdef USE_LOGGING
                Base::log.add_attribute("ID", log_id);
#endif
            };

            template<typename T>
            void set(T geometry, Identifier id) {
                Base::m_geometry = geometry;
                Base::template init<T>(geometry);
                m_id = id;
                Base::template emitSignal<reset> (Base::shared_from_this());
            };

            Identifier& getIdentifier() {
                return m_id;
            };
            void setIdentifier(Identifier id) {
                m_id = id;
#ifdef USE_LOGGING
                std::stringstream str;
                str<<id;
                log_id.set(str.str());
                BOOST_LOG(Base::log)<<"Identifyer set: "<<id;
#endif
            };
        };

        struct Geometry3D : public mpl::if_<boost::is_same<Identifier, No_Identifier>,
                detail::Geometry<Sys, Geometry3D, Typelist, GeomSignal, 3>, Geometry3D_id<Geometry3D> >::type {

            typedef typename mpl::if_<boost::is_same<Identifier, No_Identifier>,
                    detail::Geometry<Sys, Geometry3D, Typelist, GeomSignal, 3>,
                    Geometry3D_id<Geometry3D> >::type base;

            template<typename T>
            Geometry3D(T geometry, Sys& system) : base(geometry, system) { };
        };


        template<typename Derived>
        class Constraint3D_id : public detail::Constraint<Sys, Derived, ConsSignal, MES, Geometry3D> {

            typedef detail::Constraint<Sys, Derived, ConsSignal, MES, Geometry3D> base;
            Identifier m_id;
        public:
            Constraint3D_id(Sys& system, Geom f, Geom s) : base(system, f, s) {};

            Identifier& getIdentifier() {
                return m_id;
            };
            void setIdentifier(Identifier id) {
                m_id = id;
            };
        };

        struct Constraint3D : public mpl::if_<boost::is_same<Identifier, No_Identifier>,
                detail::Constraint<Sys, Constraint3D, ConsSignal, MES, Geometry3D>,
                Constraint3D_id<Constraint3D> >::type {

            typedef typename mpl::if_<boost::is_same<Identifier, No_Identifier>,
                    detail::Constraint<Sys, Constraint3D, ConsSignal, MES, Geometry3D>,
                    Constraint3D_id<Constraint3D> >::type base;

            Constraint3D(Sys& system, Geom first, Geom second) : base(system, first, second) { };
        };

        typedef mpl::vector<Geometry3D, Constraint3D> objects;

        struct inheriter_base {
            inheriter_base() {
                m_this = ((Sys*) this);
            };

            Geom drag_point, drag_goal;
            Cons drag_constraint;

            template<typename T>
            Geom createGeometry3D(T geom) {

                Geom g(new Geometry3D(geom, * ((Sys*) this)));
                fusion::vector<LocalVertex, GlobalVertex> res = m_this->m_cluster.addVertex();
                m_this->m_cluster.template setObject<Geometry3D> (fusion::at_c<0> (res), g);
                g->template setProperty<vertex_prop>(fusion::at_c<1>(res));
                m_this->template objectVector<Geometry3D>().push_back(g);
                return g;
            };

            template<typename T1>
            Cons createConstraint3D(Geom first, Geom second, T1 constraint1) {

                //build a constraint vector
                typedef mpl::vector<> cvec;
                typedef typename mpl::if_< mpl::is_sequence<T1>,
                        typename mpl::fold< T1, cvec, mpl::push_back<mpl::_1,mpl::_2> >::type,
                        mpl::vector<T1> >::type cvec1;

                //make a fusion sequence to hold the objects (as they hold the options)
                typedef typename fusion::result_of::as_vector<cvec1>::type covec;
                //set the objects
                covec cv;
                fusion::at_c<0>(cv) = constraint1;

                //now create the constraint
                Cons c(new Constraint3D(*m_this, first, second));
                //set the type and values
                c->template initialize<cvec1>(cv);

                //add it to the clustergraph
                fusion::vector<LocalEdge, GlobalEdge, bool, bool> res;
                res = m_this->m_cluster.addEdge(first->template getProperty<vertex_prop>(),
                                                second->template getProperty<vertex_prop>());
                if(!fusion::at_c<2>(res))  {
                    Cons rc;
                    return rc; //TODO: throw
                };
                m_this->m_cluster.template setObject<Constraint3D> (fusion::at_c<1> (res), c);
                //add the coresbondig edge to the constraint
                c->template setProperty<edge_prop>(fusion::at_c<1>(res));
                //store the constraint in general object vector of main system
                m_this->template objectVector<Constraint3D>().push_back(c);

                return c;
            };

        protected:
            Sys* m_this;
        };

        struct inheriter_noid : public inheriter_base {

        protected:
            using inheriter_base::m_this;

        public:
            //only point draging up to now
            bool startPointDrag(Geom g) {
                /*

                    inheriter_base::drag_point = g;
                    inheriter_base::drag_goal.reset();*/
            };

            template<typename T>
            void pointDrag(T point) {
                /*
                    BOOST_MPL_ASSERT((boost::is_same< typename geometry_traits<T>::tag, typename tag::point3D>));
                    if(!inheriter_base::drag_goal) {
                        inheriter_base::drag_goal = this->createGeometry3D(point);
                        inheriter_base::drag_constraint = this->template createConstraint3D<Distance3D>(inheriter_base::drag_point, inheriter_base::drag_goal, 0);
                    }
                    inheriter_base::drag_goal->set(point);
                    this->solve();*/
            };
            void finishPointDrag() {
                /*
                    //TODO:remove constraints and drag goal
                    inheriter_base::drag_goal.reset();
                    inheriter_base::drag_constraint.reset();*/
            };
        };

        struct inheriter_id : public inheriter_base {

        protected:
            using inheriter_base::m_this;

        public:
            template<typename T>
            Geom createGeometry3D(T geom, Identifier id) {
                Geom g = inheriter_base::createGeometry3D(geom);
                g->setIdentifier(id);
                return g;
            };

            template<typename T1>
            Cons createConstraint3D(Identifier id, Geom first, Geom second, T1 constraint1) {

                Cons c = inheriter_base::createConstraint3D(first, second, constraint1);
                c->setIdentifier(id);
                return c;
            };


            bool hasGeometry3D(Identifier id) {
                if(getGeometry3D(id)) return true;
                return false;
            };

            Geom getGeometry3D(Identifier id) {
                std::vector< Geom >& vec = inheriter_base::m_this->template objectVector<Geometry3D>();
                typedef typename std::vector<Geom>::iterator iter;
                for(iter it=vec.begin(); it!=vec.end(); it++) {
                    if(compare_traits<Identifier>::compare((*it)->getIdentifier(), id)) return *it;
                };
                return Geom();
            };

            bool hasConstraint3D(Identifier id) {
                if(getConstraint3D(id)) return true;
                return false;
            };

            Cons getConstraint3D(Identifier id) {
                std::vector< Cons >& vec = inheriter_base::m_this->template objectVector<Constraint3D>();
                typedef typename std::vector<Cons>::iterator iter;
                for(iter it=vec.begin(); it!=vec.end(); it++) {
                    if(compare_traits<Identifier>::compare((*it)->getIdentifier(), id)) return *it;
                };
                return Cons();
            };

            //only point draging up to now
            bool startPointDrag(Identifier id) {
                /*
                                inheriter_base::drag_point = getGeometry3D(id);
                                inheriter_base::drag_goal.reset();*/
            };

            template<typename T>
            void pointDrag(T point) {
                /*BOOST_MPL_ASSERT((boost::is_same< typename geometry_traits<T>::tag, typename tag::point3D>));
                if(!inheriter_base::drag_goal) {
                    inheriter_base::drag_goal = this->createGeometry3D(point, "drag_goal");
                    inheriter_base::drag_constraint = this->template createConstraint3D<Fix3D>("drag_constraint", inheriter_base::drag_point, inheriter_base::drag_goal, 0);
                }
                inheriter_base::drag_goal->set(point, "drag_goal");
                //inheriter_base::drag_goal->m_parameterCount=0;
                ((Sys*) this)->solve();*/
            };
            void finishPointDrag() {
                /*
                    //TODO:remove constraints and drag goal
                    inheriter_base::drag_goal.reset();
                    inheriter_base::drag_constraint.reset();*/
            };
        };

        struct inheriter : public mpl::if_<boost::is_same<Identifier, No_Identifier>, inheriter_noid, inheriter_id>::type {};



        struct math_prop {
            typedef cluster_property kind;
            typedef details::ClusterMath<Sys> type;
        };
        struct fix_prop {
            typedef cluster_property kind;
            typedef bool type;
        };
        struct vertex_prop {
            typedef Geometry3D kind;
            typedef GlobalVertex type;
        };
        struct edge_prop {
            typedef Constraint3D kind;
            typedef GlobalEdge type;
        };

        typedef mpl::vector<vertex_prop, edge_prop, math_prop, fix_prop>  properties;

        static void system_init(Sys& sys) {
            sys.m_sheduler.addProcessJob(new SystemSolver());
        };
    };
};

namespace details {
//allow direct access to the stored geometry in a Geometry3D, copyed from boost variant get
template <typename T>
struct get_visitor {
private:

    typedef typename boost::add_pointer<T>::type pointer;
    typedef typename boost::add_reference<T>::type reference;

public:

    typedef pointer result_type;

public:
    pointer operator()(reference operand) const   {
        return boost::addressof(operand);
    }

    template <typename U>
    pointer operator()(const U&) const  {
        return static_cast<pointer>(0);
    }
};
}

template<typename T, typename G>
typename boost::add_reference<T>::type get(G geom) {

    typedef typename boost::add_pointer<T>::type T_ptr;
    details::get_visitor<T> v;
    T_ptr result = geom->apply(v);

    //if (!result)
    //TODO:throw bad_get();
    return *result;
};

}//dcm

#endif //GCM_GEOMETRY3D_H







