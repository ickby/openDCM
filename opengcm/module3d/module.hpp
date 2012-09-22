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

#ifndef GCM_MODULE_3D_H
#define GCM_MODULE_3D_H

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

#include "opengcm/Core"
#include "opengcm/core/object.hpp"
#include "opengcm/core/clustergraph.hpp"
#include "opengcm/core/sheduler.hpp"
#include "opengcm/core/traits.hpp"
#include "opengcm/core/geometry.hpp"
#include "geometry.hpp"
#include "coincident.hpp"
#include "distance.hpp"
#include "parallel.hpp"
#include "angle.hpp"
#include "dof.hpp"
#include "fix.hpp"


#define maxfak 1.2
#define minfak 0.8

static int counter = 0;

namespace mpl = boost::mpl;

namespace gcm {

namespace details {

enum { cluster3D = 100};

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
    typename Kernel::Vector3	m_original_translation;
    typename Kernel::Vector3Map	m_normQ;
    typename Kernel::Vector3 	m_origNormQ;

    int m_rot_offset, m_trans_offset;
    int count;
    bool init;

public:
    typename Kernel::Vector3Map m_translation;
    //shift scale stuff
    Scalar xmin, xmax, ymin, ymax, zmin, zmax;
    typename Kernel::Vector3 midpoint, m_shift;
    Scalar m_scale;

public:
    ClusterMath() : m_normQ(NULL), m_translation(NULL), init(false) {

        m_quaternion = typename Kernel::Quaternion(1,2,3,4);
        m_quaternion.normalize();
        m_shift.setZero();
        count = counter;
        counter++;
        m_scale = 1.;
        xmin=1e10;
        xmax=-1e10;
        ymin=1e10;
        ymax=-1e10;
        zmin=1e10;
        zmax=-1e10;
    };

    void setParameterOffset(int roff, int toff) {
        m_rot_offset = roff;
        m_trans_offset = toff;
    };
    int getRotationOffset() {
        return m_rot_offset;
    };
    int getTranslationOffset() {
        return m_trans_offset;
    };

    void setRotationMap(typename Kernel::Matrix3Map& map, typename Kernel::Matrix39Map& diffmap) {
        new(&map) typename Kernel::Matrix3Map(&m_rotation(0,0),3,3);
        new(&diffmap) typename Kernel::Matrix39Map(&m_diffrot(0,0));
    };
    void setTranslationMap(typename Kernel::Vector3Map& map) {
        new(&map) typename Kernel::Vector3Map(&m_translation(0));
    };
    void setShiftMap(typename Kernel::Vector3Map& map) {
        new(&map) typename Kernel::Vector3Map(&m_shift(0));
    };
    typename Kernel::Vector3Map& getNormQuaternionMap() {
        return m_normQ;
    };
    typename Kernel::Vector3Map& getTranslationMap() {
        return m_translation;
    };
    void initMaps() {
        const Scalar s = std::acos(m_quaternion.w())/std::sin(std::acos(m_quaternion.w()));
        m_normQ = m_quaternion.vec()*s;
        m_translation = m_original_translation + m_quaternion.toRotationMatrix()*m_shift;
	init = true;
	m_scale = 1.;
        xmin=1e10;
        xmax=-1e10;
        ymin=1e10;
        ymax=-1e10;
        zmin=1e10;
        zmax=-1e10;
	midpoint.setZero();
	m_shift.setZero();
    };

    typename Kernel::Quaternion& getQuaternion() {
        return m_quaternion;
    };
    typename Kernel::Vector3& getTranslation() {
        return m_original_translation;
    };
    void setShift(typename Kernel::Vector3 s) {
        m_shift = s;
        //we remove shift from the local geometries, therefore we have to add it here
        //to not change the global position
        if(init) m_translation += m_quaternion.toRotationMatrix()*m_shift;
    };
    void setScale(Scalar s) {
        m_scale = s;
        if(init)m_translation *= s;
    };

    void finishCalculation() {
        const Scalar norm = m_normQ.norm();
        const Scalar fac = std::sin(norm)/norm;
        m_quaternion = typename Kernel::Quaternion(std::cos(norm), m_normQ(0)*fac, m_normQ(1)*fac, m_normQ(2)*fac);
        m_quaternion.normalize();
        m_original_translation = m_translation/m_scale - m_quaternion.toRotationMatrix()*m_shift;

        //needed to allow a correct global calculation in cluster geometries after this finish
        m_shift.setZero();
        m_translation = m_original_translation;
	
	init=false;
    };

    void recalculate() {

        //get the Quaternion for the norm quaternion form and calculate the rotation matrix
        const Scalar norm = m_normQ.norm();
        const Scalar fac = std::sin(norm)/norm;

        typename Kernel::Quaternion Q(std::cos(norm), m_normQ(0)*fac, m_normQ(1)*fac, m_normQ(2)*fac);
        Q.normalize(); //not needed, just to avoid rounding errors
        if(Kernel::isSame(norm, 0)) {
            Q.setIdentity();
            m_rotation.setIdentity();
            m_diffrot.setZero();
            return;
        };
        m_rotation = Q.toRotationMatrix();

        /* now calculate the gradient quaternions and calculate the diff rotation matrices
         * m_normQ = (a,b,c)
         * n = ||m_normQ||
         *
         * Q = (a/n sin(n), b/n sin(n), c/n sin(n), cos(n))
         */

        //n=||m_normQ||, sn = sin(n)/n, sn3 = sin(n)/n^3, cn = cos(n)/n, divn = 1/n;
        const Scalar n    = m_normQ.norm();
        const Scalar sn   = std::sin(n)/n;
        const Scalar mul  = (std::cos(n)-sn)/std::pow(n,2);

        //dxa = dx/da
        const Scalar dxa = sn + std::pow(m_normQ(0),2)*mul;
        const Scalar dxb = m_normQ(0)*m_normQ(1)*mul;
        const Scalar dxc = m_normQ(0)*m_normQ(2)*mul;

        const Scalar dya = m_normQ(1)*m_normQ(0)*mul;
        const Scalar dyb = sn + std::pow(m_normQ(1),2)*mul;
        const Scalar dyc = m_normQ(1)*m_normQ(2)*mul;

        const Scalar dza = m_normQ(2)*m_normQ(0)*mul;
        const Scalar dzb = m_normQ(2)*m_normQ(1)*mul;
        const Scalar dzc = sn + std::pow(m_normQ(2),2)*mul;

        const Scalar dwa = -sn*m_normQ(0);
        const Scalar dwb = -sn*m_normQ(1);
        const Scalar dwc = -sn*m_normQ(2);

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

struct m3d {}; 	//base of module3d::type to allow other modules check for it

template<typename Typelist, typename Identifier = No_Identifier>
struct Module3D {

    template<typename Sys>
    struct type : m3d {
        class Constraint3D;
        class Geometry3D;
        typedef boost::shared_ptr<Geometry3D> Geom;
        typedef boost::shared_ptr<Constraint3D> Cons;

        typedef mpl::map< mpl::pair<reset, boost::function<void (Geom) > > >  GeomSignal;

        struct MES  : public system_traits<Sys>::Kernel::MappedEquationSystem {

            typedef typename system_traits<Sys>::Cluster Cluster;
            Cluster& m_cluster;

            MES(Cluster& cl, int par, int rpar, int tpar, int eqn)
                : system_traits<Sys>::Kernel::MappedEquationSystem(par, rpar, tpar, eqn),
                  m_cluster(cl) {};

            virtual void recalculate() {

                //first calculate all clusters
                typedef typename Cluster::cluster_iterator citer;
                std::pair<citer, citer> cit = m_cluster.clusters();
                for(; cit.first != cit.second; cit.first++) {

                    if(!(*cit.first).second->template getClusterProperty<fix_prop>()) {
                        (*cit.first).second->template getClusterProperty<math_prop>().recalculate();

                        //now with the new rotation matrix we calculate all geometries in that cluster
                        std::vector<Geom>& vec = (*cit.first).second->template getClusterProperty<gmap_prop>();
                        typedef typename std::vector<Geom>::iterator iter;

                        for(iter it = vec.begin(); it != vec.end(); it++)
                            (*it)->recalculate(system_traits<Sys>::Kernel::MappedEquationSystem::Scaling);
                    }

                };
                //TODO:Scale parameters outside cluster

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
                            (*oit.first)->calculate();
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
                //std::cout<<"hmmm"<<std::endl;
                for(; cit.first != cit.second; cit.first++) {

                    if((*cit.first).second->template getClusterProperty<changed_prop>() &&
                            (*cit.first).second->template getClusterProperty<type_prop>() == details::cluster3D)
                        solveCluster(*(*cit.first).second, sys);
                }

                int params=0, trans_params=0, rot_params=0, constraints=0;
                typename Kernel::number_type scale = 1;

                //get the ammount of parameters and constraint equations we need
                typedef typename boost::graph_traits<Cluster>::vertex_iterator iter;
                std::pair<iter, iter>  it = boost::vertices(cluster);
                for(; it.first != it.second; it.first++) {

                    //when cluster and not fixed it has trans and rot parameter
                    if(cluster.isCluster(*it.first)) {
                        if(!cluster.template getSubclusterProperty<fix_prop>(*it.first)) {
                            trans_params += 3;
                            rot_params += 3;
                        }
                    } else {
                        params += cluster.template getObject<Geometry3D>(*it.first)->m_parameterCount;
                    };
                }

                typedef typename boost::graph_traits<Cluster>::edge_iterator e_iter;
                std::pair<e_iter, e_iter>  e_it = boost::edges(cluster);
                for(; e_it.first != e_it.second; e_it.first++)
                    constraints += cluster.getGlobalEdgeCount(*e_it.first);

                //initialise the system with now known size
                MES mes(cluster, params, rot_params, trans_params, constraints);

                //iterate all geometrys again and set the needed maps
                it = boost::vertices(cluster);
                for(; it.first != it.second; it.first++) {

                    if(cluster.isCluster(*it.first)) {
                        Cluster& c = cluster.getVertexCluster(*it.first);
                        details::ClusterMath<Sys>& cm =  c.template getClusterProperty<math_prop>();
                        //only get maps and propagate downstream if not fixed
                        if(!c.template getClusterProperty<fix_prop>()) {
                            //set norm Quaternion as map to the parameter vector
                            int offset = mes.setParameterMap(Rotation, cm.getNormQuaternionMap());
                            //set translation as map to the parameter vector
                            int transoffset = mes.setParameterMap(Translation, cm.getTranslationMap());
                            //write initail values to the parameter maps
                            //remember the parameter offset as all downstream geometry must use this offset
                            cm.setParameterOffset(offset, transoffset);
                            //wirte initial values
                            cm.initMaps();
                        }

                        //map all geometrie within that cluster to it's rotation matrix
                        //for collecting all geometries which need updates
                        std::vector<Geom>& vec = c.template getClusterProperty<gmap_prop>();
                        vec.clear();
                        //to allow a corect calculation of geometries toplocal value we need the quaternion
                        //which transforms from toplevel to this "to be solved" cluster and the aquivalent translation
                        typename Kernel::Quaternion q(1,0,0,0);
                        typename Kernel::Vector3 t(0,0,0);
                        mapClusterDownstreamGeometry(c, cm, vec, q, t);


                    } else {
                        Geom g = cluster.template getObject<Geometry3D>(*it.first);
                        int offset = mes.setParameterMap(Anything, g->m_parameterCount, g->getParameterMap());
                        g->m_parameter_offset = offset;
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
                        if(c) {

                            //when in cluster, there are 6 clusterparameter we differentiat for, if not we differentiat
                            //for every parameter in the geometry;
                            int equation = mes.setResidualMap(c->m_residual);
                            if(c->first->getClusterMode()) {
                                if(!c->first->isClusterFixed()) {
                                    mes.setJacobiMap(equation, c->first->m_trans_offset, 3, c->m_trans_diff_first);
                                    mes.setJacobiMap(equation, c->first->m_rot_offset, 3, c->m_rot_diff_first);
                                }
                            } else  mes.setJacobiMap(equation, c->first->m_parameter_offset, c->first->m_parameterCount, c->m_diff_first);

                            if(c->second->getClusterMode()) {
                                if(!c->second->isClusterFixed()) {
                                    mes.setJacobiMap(equation, c->second->m_trans_offset, 3, c->m_trans_diff_second);
                                    mes.setJacobiMap(equation, c->second->m_rot_offset, 3, c->m_rot_diff_second);
                                }
                            } else mes.setJacobiMap(equation, c->second->m_parameter_offset, c->first->m_parameterCount, c->m_diff_second);

                            //mes.Scaling(equation) = c->getEquationScaling();
                        }
                        //TODO: else throw (as every global edge was counted as one equation)
                    }
                }

                //get the maximal scale
                Scalar sc = 0;
                for(cit = cluster.clusters(); cit.first != cit.second; cit.first++) {
                    const Scalar s = scaleShiftCluster(*(*cit.first).second);
                    sc = (s>sc) ? s : sc;
                }
                //Base::Console().Message("Scale is %f\n", sc);
                if(!Kernel::isSame(sc,0)) {
                    for(cit = cluster.clusters(); cit.first != cit.second; cit.first++)
                        applyShiftCluster(*(*cit.first).second, sc);
                    mes.Scaling = maxfak/sc;
                }
                //scaling needs to be 1 (all shifts are at all theri points)
                else {
                    for(cit = cluster.clusters(); cit.first != cit.second; cit.first++)
                        applyShiftCluster(*(*cit.first).second, 1.);
                    mes.Scaling = 1.;
                }

                //now it's time to solve
                Kernel::solve(mes);

                std::cout<<"Residual after solving: "<<mes.Residual.norm()<<std::endl;

                //now go to all relevant geometries and clusters and write the values back
                it = boost::vertices(cluster);
                for(; it.first != it.second; it.first++) {

                    if(cluster.isCluster(*it.first)) {
                        if(!cluster.template getSubclusterProperty<fix_prop>(*it.first)) {
                            Cluster& c = cluster.getVertexCluster(*it.first);

                            c.template getClusterProperty<math_prop>().finishCalculation();
                            std::vector<Geom>& vec = c.template getClusterProperty<gmap_prop>();
                            for(typename std::vector<Geom>::iterator vit = vec.begin(); vit != vec.end(); vit++)
                                (*vit)->finishCalculation();
                        }
                    } else cluster.template getObject<Geometry3D>(*it.first)->finishCalculation();
                }

                //we have solved this cluster
                cluster.template setClusterProperty<changed_prop>(false);

            };


            void mapClusterDownstreamGeometry(Cluster& cluster,
                                              details::ClusterMath<Sys>& cm,
                                              std::vector<Geom>& vec,
                                              typename Kernel::Quaternion& q,
                                              typename Kernel::Vector3& t) {
                //all geometry within that cluster needs to be mapped to the provided rotation matrix (in cm)
                //also the geometries toplocal value needs to be set so that it matches this cm
                typename Kernel::Quaternion nq = q*cluster.template getClusterProperty<math_prop>().getQuaternion();
                typename Kernel::Vector3 nt = t+cluster.template getClusterProperty<math_prop>().getTranslation();

                //get all vertices and map the geometries if existend
                typedef typename boost::graph_traits<Cluster>::vertex_iterator iter;
                std::pair<iter, iter>  it = boost::vertices(cluster);
                for(; it.first != it.second; it.first++) {
                    Geom g = cluster.template getObject<Geometry3D>(*it.first);
                    if(g) {
                        if(!cluster.template getClusterProperty<fix_prop>()) {
                            //allow iteration over all maped geometries
                            vec.push_back(g);
                            //map rotation and diffrotation from cluster to geometry
                            cm.setRotationMap(g->getRotationMap(), g->getDiffRotationMap());
                            //map translation from cluster to geometry
                            cm.setTranslationMap(g->getTranslationMap());
                            //map shift from cluster to geometry
                            cm.setShiftMap(g->getShiftMap());
                            //set the offsets so that geometry knows where it is in the parameter map
                            g->m_rot_offset = cm.getRotationOffset();
                            g->m_trans_offset = cm.getTranslationOffset();
                        }
                        //calculate the appropriate local values
                        g->transformInverse(nq.conjugate().toRotationMatrix(), -nt);

                        //position and offset of the parameters must be set to the clusters values
                        g->setClusterMode(true, cluster.template getClusterProperty<fix_prop>());
                    }
                }

                //go downstream and map
                typedef typename Cluster::cluster_iterator citer;
                std::pair<citer, citer> cit = cluster.clusters();
                for(; cit.first != cit.second; cit.first++)
                    mapClusterDownstreamGeometry(*(*cit.first).second, cm, vec, nq, nt);
                //TODO: if one subcluster is fixed the hole cluster should be too, as there are no
                //	dof's remaining between parts and so nothing can be moved when one part is fixed.

            };

            Scalar scaleShiftCluster(Cluster& cluster) {

                //first get the bonding box to get the center pof points
                std::vector<Geom>& vec = cluster.template getClusterProperty<gmap_prop>();
                details::ClusterMath<Sys>& math = cluster.template getClusterProperty<math_prop>();

                std::stringstream str;
                str<<"Start ScaleShift"<<std::endl;
                //only one geometry scale = norm
                if(vec.empty()) return 0.; //should never happen...
                if(vec.size() == 1) {
                    str<<"single geometry part: "<<vec[0]->getBigPoint().transpose()<<std::endl;
                    //Base::Console().Message("%s",str.str().c_str());
                    return vec[0]->getBigPoint().norm();
                }

                typedef typename std::vector<Geom>::iterator iter;
                for(iter it = vec.begin(); it != vec.end(); it++) {
                    typename Kernel::Vector3 v = (*it)->getBigPoint();
                    math.xmin = (v(0)<math.xmin) ? v(0) : math.xmin;
                    math.xmax = (v(0)<math.xmax) ? math.xmax : v(0);
                    math.ymin = (v(1)<math.ymin) ? v(1) : math.ymin;
                    math.ymax = (v(1)<math.ymax) ? math.ymax : v(1);
                    math.zmin = (v(2)<math.zmin) ? v(2) : math.zmin;
                    math.zmax = (v(2)<math.zmax) ? math.zmax : v(2);
                    str<<"Testet point: "<<v.transpose()<<std::endl;
                };

                //now calculate the midpoint and use it as shift
                math.midpoint << math.xmin+math.xmax, math.ymin+math.ymax, math.zmin+math.zmax;
                math.midpoint /= 2.;
                str<<"Midpoint:"<<math.midpoint.transpose()<<std::endl;
                //the bounding box corner is the max allowed distance
                typename Kernel::Vector3 max(math.xmax, math.ymax, math.zmax);
                Scalar maxscale = (max-math.midpoint).norm();

                //the maxscale is ||point||*scale=1.5, all other points are allowed to be in the
                //range to ||point||*scale=0.5. therefore ||point|| > maxscale/3 must be ensured
                bool inFrame = true;
                for(iter it = vec.begin(); (it != vec.end()) && inFrame; it++)
                    inFrame = (((*it)->getBigPoint()-math.midpoint).norm() > maxscale*minfak/maxfak);

                if(!inFrame) str<<"Point not in right shift range"<<std::endl;
                //all points are in frame, we are done here
                //if(inFrame) return 2.*maxscale/3.;
                math.m_scale =  maxscale;

                str<<"Scale: "<<math.m_scale<<std::endl;

                //set the calulated shift
                math.m_shift = math.midpoint;

                str<<"shift: "<<math.m_shift.transpose()<<std::endl<<std::endl;
                //Base::Console().Message("%s",str.str().c_str());
                return maxscale;

                //some points are to close to the origin, lets shift again
                //first get the flat boundingbox side

            };

            void applyShiftCluster(Cluster& cluster, Scalar scale) {

                details::ClusterMath<Sys>& math = cluster.template getClusterProperty<math_prop>();

                //same scaling for us all
                //math.m_scale = scale;

                //if only one point exists we extend the origin-point-line to match the scale
                std::vector<Geom>& vec = cluster.template getClusterProperty<gmap_prop>();
                if(vec.size()==1) {
                    typename Kernel::Vector3 v = vec[0]->getBigPoint();
                    const Scalar fak = 1. - scale/v.norm();
                    if(Kernel::isSame(v.norm(),0))
                        math.m_shift << scale, 0, 0;
                    else math.m_shift = fak*v;

                    math.setShift(math.m_shift);
                    math.setScale(maxfak/scale);

                    std::stringstream str;
                    str<<"single shift: "<<math.m_shift.transpose()<<std::endl<<std::endl;
                    //Base::Console().Message("%s",str.str().c_str());
                    return;
                };


                //if this is our scale then just applie the midpoint as shift (already done)
                if(Kernel::isSame(scale, math.m_scale)) {
                    math.setShift(math.m_shift);
                    math.setScale(maxfak/scale);
                    return;
                }

                //now it gets more complicated. bahh. the most outer point should be at
                //distance scale. of course we need the smallest heigh at midpoint to add
                //the offset to
                Scalar xh = math.xmax-math.xmin;
                Scalar yh = math.ymax-math.ymin;
                Scalar zh = math.zmax-math.zmin;

                if((xh<=yh) && (xh<=zh)) {
                    const Scalar a = std::sqrt(std::pow(math.ymax-math.midpoint(1),2)
                                               +  std::pow(math.zmax-math.midpoint(2),2));
                    const Scalar b = std::sqrt(std::pow(scale,2) - std::pow(a,2));

                    //applie extra heigh to midpoint
                    math.m_shift(0) += b-xh/2;
                } else if((yh<xh) && (yh<zh)) {
                    const Scalar a = std::sqrt(std::pow(math.xmax-math.midpoint(0),2)
                                               +  std::pow(math.zmax-math.midpoint(2),2));
                    const Scalar b = std::sqrt(std::pow(scale,2) - std::pow(a,2));

                    //applie extra heigh to midpoint
                    math.m_shift(1) += b-yh/2;
                } else {
                    const Scalar a = std::sqrt(std::pow(math.xmax-math.midpoint(0),2)
                                               +  std::pow(math.ymax-math.midpoint(1),2));
                    const Scalar b = std::sqrt(std::pow(scale,2) - std::pow(a,2));

                    //applie extra heigh to midpoint
                    math.m_shift(2) += b-zh/2;
                }

                math.setShift(math.m_shift);
                math.setScale(maxfak/scale);

                std::stringstream str;
                str<<"changed shift: "<<math.m_shift.transpose()<<std::endl<<std::endl;
                //Base::Console().Message("%s",str.str().c_str());

            };

        };

	template<typename Derived>
        class Geometry3D_id : public detail::Geometry<Sys, Derived, Typelist, GeomSignal> {

	    typedef detail::Geometry<Sys, Derived, Typelist, GeomSignal> Base;
	    
            Identifier m_id;
        public:
            template<typename T>
            Geometry3D_id(T geometry, Sys& system) : Base(geometry, system) { };

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
            };
        };

        struct Geometry3D : public mpl::if_<boost::is_same<Identifier, No_Identifier>,
                detail::Geometry<Sys, Geometry3D, Typelist, GeomSignal>, Geometry3D_id<Geometry3D> >::type {

            typedef typename mpl::if_<boost::is_same<Identifier, No_Identifier>,
                    detail::Geometry<Sys, Geometry3D, Typelist, GeomSignal>,
		    Geometry3D_id<Geometry3D> >::type base;

            template<typename T>
            Geometry3D(T geometry, Sys& system) : base(geometry, system) { };
        };


        //type erasure container for constraints
        class Constraint3D_base : public Object<Sys, Constraint3D, mpl::map<> > {

            typedef typename system_traits<Sys>::Kernel Kernel;
            typedef typename Kernel::number_type Scalar;
            typedef typename Kernel::DynStride DS;

        public:
            Constraint3D_base(Sys& system, Geom f, Geom s) : Object<Sys, Constraint3D, mpl::map<> > (system),
                first(f), second(s), content(0), m_diff_first(NULL,0,DS(0,0)), m_diff_second(NULL,0,DS(0,0)),
                m_rot_diff_first(NULL,0,DS(0,0)), m_rot_diff_second(NULL,0,DS(0,0)),
                m_trans_diff_first(NULL,0,DS(0,0)), m_trans_diff_second(NULL,0,DS(0,0)), m_residual(NULL,0,DS(0,0))	{

                cf = first->template connectSignal<reset> (boost::bind(&Constraint3D::geometryReset, this, _1));
                cs = second->template connectSignal<reset> (boost::bind(&Constraint3D::geometryReset, this, _1));
            };

            ~Constraint3D_base()  {
                delete content;
                first->template disconnectSignal<reset>(cf);
                second->template disconnectSignal<reset>(cs);
            }

        protected:
            template< typename creator_type>
            void resetType(creator_type& c) {
                boost::apply_visitor(c, first->m_geometry, second->m_geometry);
                content = c.p;
                if(c.need_swap) first.swap(second);
            };

            Scalar getEquationScaling() {
                return content->getEquationScaling(first->m_toplocal, second->m_toplocal);
            };

            Scalar calculate() {

                m_residual(0) = content->calculate(first->m_parameter, second->m_parameter);

                //now see which way we should calculate the gradient (may be diffrent for both geometries)
                if(first->m_parameterCount) {
                    if(first->getClusterMode()) {
                        //when the cluster is fixed no maps are set as no parameters exist.
                        if(!first->isClusterFixed()) {

                            //cluster mode, so we do a full calculation with all 3 rotation diffparam vectors
                            for(int i=0; i<3; i++) {
                                typename Kernel::VectorMap block(&first->m_diffparam(0,i),first->m_parameterCount,1, DS(1,1));
                                m_rot_diff_first(i) = content->calculateGradientFirst(first->m_parameter,
                                                      second->m_parameter, block);
                            }
                            //now the translation stuff
                            for(int i=3; i<6; i++) {
                                typename Kernel::VectorMap block(&first->m_diffparam(0,i),first->m_parameterCount,1, DS(1,1));
                                m_trans_diff_first(i-3) = content->calculateGradientFirst(first->m_parameter,
                                                          second->m_parameter, block);
                            }
                        }
                    } else {
                        //not in cluster, so allow the constraint to optimize the gradient calculation
                        content->calculateGradientFirstComplete(first->m_parameter, second->m_parameter, m_diff_first);
                    }
                }
                if(second->m_parameterCount) {
                    if(second->getClusterMode()) {
                        if(!second->isClusterFixed()) {

                            //cluster mode, so we do a full calculation with all 3 rotation diffparam vectors
                            for(int i=0; i<3; i++) {
                                typename Kernel::VectorMap block(&second->m_diffparam(0,i),second->m_parameterCount,1, DS(1,1));
                                m_rot_diff_second(i) = content->calculateGradientSecond(first->m_parameter,
                                                       second->m_parameter, block);
                            }
                            //translational gradients
                            for(int i=3; i<6; i++) {
                                typename Kernel::VectorMap block(&second->m_diffparam(0,i),second->m_parameterCount,1, DS(1,1));
                                m_trans_diff_second(i-3) = content->calculateGradientSecond(first->m_parameter,
                                                           second->m_parameter, block);
                            }
                        }
                    } else {
                        //not in cluster, so allow the constraint to optimize the gradient calculation
                        content->calculateGradientSecondComplete(first->m_parameter, second->m_parameter, m_diff_second);
                    }

                }
            };

            void geometryReset(Geom g) {
                placeholder* p = content->resetConstraint(first, second);
                delete content;
                content = p;
            };

            struct placeholder  {

                virtual ~placeholder() {}
                virtual placeholder* resetConstraint(Geom first, Geom second) const = 0;

                virtual Scalar calculate(typename Kernel::VectorMap&, typename Kernel::VectorMap&) = 0;
                virtual Scalar calculateGradientFirst(typename Kernel::VectorMap& param1,
                                                      typename Kernel::VectorMap& param2,
                                                      typename Kernel::VectorMap& dparam1) = 0;
                virtual Scalar calculateGradientSecond(typename Kernel::VectorMap& param1,
                                                       typename Kernel::VectorMap& param2,
                                                       typename Kernel::VectorMap& dparam2) = 0;
                virtual void calculateGradientFirstComplete(typename Kernel::VectorMap& param1,
                        typename Kernel::VectorMap& param2,
                        typename Kernel::VectorMap& grad) = 0;
                virtual void calculateGradientSecondComplete(typename Kernel::VectorMap& param1,
                        typename Kernel::VectorMap& param2,
                        typename Kernel::VectorMap& grad) = 0;
                virtual Scalar getEquationScaling(typename Kernel::Vector& local1,
                                                  typename Kernel::Vector& local2) = 0;
            };

            template< template<typename, typename, typename> class T1, typename T2, typename T3>
            struct holder : public placeholder  {

                holder(const T1<Kernel, T2, T3> & value)
                    : held(value)   {}

                virtual Scalar calculate(typename Kernel::VectorMap& f, typename Kernel::VectorMap& s) {
                    return held.calculate(f,s);
                };
                virtual Scalar calculateGradientFirst(typename Kernel::VectorMap& param1,
                                                      typename Kernel::VectorMap& param2,
                                                      typename Kernel::VectorMap& dparam1) {
                    return held.calculateGradientFirst(param1, param2, dparam1);
                };
                virtual Scalar calculateGradientSecond(typename Kernel::VectorMap& param1,
                                                       typename Kernel::VectorMap& param2,
                                                       typename Kernel::VectorMap& dparam2) {
                    return held.calculateGradientSecond(param1, param2, dparam2);
                };
                virtual void calculateGradientFirstComplete(typename Kernel::VectorMap& param1,
                        typename Kernel::VectorMap& param2,
                        typename Kernel::VectorMap& grad) {
                    held.calculateGradientFirstComplete(param1, param2, grad);
                };
                virtual void calculateGradientSecondComplete(typename Kernel::VectorMap& param1,
                        typename Kernel::VectorMap& param2,
                        typename Kernel::VectorMap& grad) {
                    held.calculateGradientSecondComplete(param1, param2, grad);
                };
                virtual Scalar getEquationScaling(typename Kernel::Vector& local1,
                                                  typename Kernel::Vector& local2) {
                    return held.getEquationScaling(local1, local2);
                };

                virtual placeholder* resetConstraint(Geom first, Geom second) const {
                    creator<T1> creator;
                    boost::apply_visitor(creator, first->m_geometry, second->m_geometry);
                    if(creator.need_swap) first.swap(second);
                    return creator.p;
                };

                T1<Kernel, T2, T3>  held;
            };

            template< template<typename,typename,typename> class T >
            struct creator : public boost::static_visitor<void> {

                template<typename T1, typename T2>
                void operator()(const T1&, const T2&) {
                    typedef tag_order< typename geometry_traits<T1>::tag, typename geometry_traits<T2>::tag > order;
                    typedef T<Kernel, typename order::first_tag, typename order::second_tag > type;
                    p = new holder< T, typename order::first_tag, typename order::second_tag > (type());
                    need_swap = order::swapt::value;
                };
                placeholder* p;
                bool need_swap;
            };

            template< template<typename,typename,typename> class T, typename Value1 >
            struct creator_1v : public boost::static_visitor<void> {

                creator_1v(Value1 v) : m_value(v) {};

                template<typename T1, typename T2>
                void operator()(const T1&, const T2&) {
                    typedef tag_order< typename geometry_traits<T1>::tag, typename geometry_traits<T2>::tag > order;
                    typedef T<Kernel, typename order::first_tag, typename order::second_tag > type;
                    p = new holder< T, typename order::first_tag, typename order::second_tag > (type(m_value));
                    need_swap = order::swapt::value;
                };
                placeholder* p;
                bool need_swap;
                Value1 m_value;
            };

            template< template<typename,typename,typename> class T, typename Value1, typename Value2 >
            struct creator_2v : public boost::static_visitor<void> {

                creator_2v(Value1 v1, Value2 v2) : m_value1(v1), m_value2(v2) {};

                template<typename T1, typename T2>
                void operator()(const T1&, const T2&) {
                    typedef tag_order< typename geometry_traits<T1>::tag, typename geometry_traits<T2>::tag > order;
                    typedef T<Kernel, typename order::first_tag, typename order::second_tag > type;
                    p = new holder< T, typename order::first_tag, typename order::second_tag > (type(m_value1, m_value2));
                    need_swap = order::swapt::value;
                };
                placeholder* p;
                bool need_swap;
                Value1 m_value1;
                Value2 m_value2;
            };


            placeholder* content;
            Geom first, second;
            Connection cf, cs;

            typename Kernel::VectorMap m_rot_diff_first, m_trans_diff_first, m_diff_first; //first geometry diff
            typename Kernel::VectorMap m_rot_diff_second, m_trans_diff_second, m_diff_second; //second geometry diff
            typename Kernel::VectorMap m_residual;

            friend class SystemSolver;
            friend class MES;

        };

        class Constraint3D_noid : public Constraint3D_base {

            using Constraint3D_base::resetType;
        public:
            Constraint3D_noid(Sys& system, Geom f, Geom s) : Constraint3D_base(system, f, s) {};

            template< template<typename,typename,typename> class T>
            void set() {
                typename Constraint3D_base::template creator<T> c;
                resetType(c);
            };
            template< template<typename,typename,typename> class T, typename Value1>
            void set(Value1 v) {
                typename Constraint3D_base::template creator_1v<T, Value1> c(v);
                resetType(c);
            };
            template< template<typename,typename,typename> class T, typename Value1, typename Value2>
            void set(Value1 v1, Value2 v2) {
                typename Constraint3D_base::template creator_2v<T, Value1, Value2> c(v1,v2);
                resetType(c);
            };
        };

        class Constraint3D_id : public Constraint3D_base {

            using Constraint3D_base::resetType;
            Identifier m_id;
        public:
            Constraint3D_id(Sys& system, Geom f, Geom s) : Constraint3D_base(system, f, s) {};

            Identifier& getIdentifier() {
                return m_id;
            };
            void setIdentitfier(Identifier id) {
                m_id = id;
            };

            template< template<typename,typename,typename> class T>
            void set(Identifier id) {
                typename Constraint3D_base::template creator<T> c;
                resetType(c);
                m_id = id;
            };
            template< template<typename,typename,typename> class T, typename Value1>
            void set(Identifier id, Value1 v) {
                typename Constraint3D_base::template creator_1v<T, Value1> c(v);
                resetType(c);
                m_id = id;
            };
            template< template<typename,typename,typename> class T, typename Value1, typename Value2>
            void set(Identifier id, Value1 v1, Value2 v2) {
                typename Constraint3D_base::template creator_2v<T, Value1, Value2> c(v1,v2);
                resetType(c);
                m_id = id;
            };
        };

        struct Constraint3D : public mpl::if_<boost::is_same<Identifier, No_Identifier>,
                Constraint3D_noid, Constraint3D_id>::type {

            typedef typename mpl::if_<boost::is_same<Identifier, No_Identifier>,
                    Constraint3D_noid, Constraint3D_id>::type base;

            Constraint3D(Sys& system, Geom f, Geom s) : base(system, f, s) { };
        };

        typedef mpl::vector<Geometry3D, Constraint3D> objects;

        struct inheriter_base {
            inheriter_base() {
                m_this = ((Sys*) this);
            };

            Geom drag_point, drag_goal;
            Cons drag_constraint;

        protected:
            Sys* m_this;

            void process_constraint(Cons c, Geom first, Geom second) {

                fusion::vector<LocalEdge, GlobalEdge, bool, bool> res;
                res = m_this->m_cluster.addEdge(first->template getProperty<vertex_prop>(),
                                                second->template getProperty<vertex_prop>());
                if(!fusion::at_c<2>(res))  {
                    return; //TODO: throw
                };
                m_this->m_cluster.template setObject<Constraint3D> (fusion::at_c<1> (res), c);
                c->template setProperty<edge_prop>(fusion::at_c<1>(res));
                m_this->template objectVector<Constraint3D>().push_back(c);
            }
        };

        struct inheriter_noid : public inheriter_base {

        protected:
            using inheriter_base::m_this;
            using inheriter_base::process_constraint;

        public:
            template<typename T>
            Geom createGeometry3D(T geom) {

                Geom g(new Geometry3D(geom, * ((Sys*) this)));
                fusion::vector<LocalVertex, GlobalVertex> res = m_this->m_cluster.addVertex();
                m_this->m_cluster.template setObject<Geometry3D> (fusion::at_c<0> (res), g);
                g->template setProperty<vertex_prop>(fusion::at_c<1>(res));
                m_this->template objectVector<Geometry3D>().push_back(g);
                return g;
            };

            template< template<typename,typename,typename> class T >
            Cons createConstraint3D(Geom first, Geom second) {

                Cons c(new Constraint3D(* ((Sys*) this), first, second));
                c->template set<T>();
                process_constraint(c, first, second);
                return c;
            };
            template< template<typename,typename,typename> class T, typename Value1 >
            Cons createConstraint3D(Geom first, Geom second, Value1 v1) {

                Cons c(new Constraint3D(* ((Sys*) this), first, second));
                c->template set<T>(v1);
                process_constraint(c, first, second);
                return c;
            };
            template< template<typename,typename,typename> class T, typename Value1, typename Value2 >
            Cons createConstraint3D(Geom first, Geom second, Value1 v1, Value2 v2) {

                Cons c(new Constraint3D(* ((Sys*) this), first, second));
                c->template set<T>(v1,v2);
                process_constraint(c, first, second);
                return c;
            };

            //only point draging up to now
            bool startPointDrag(Geom g) {

                inheriter_base::drag_point = g;
                inheriter_base::drag_goal.reset();
            };

            template<typename T>
            void pointDrag(T point) {
                BOOST_MPL_ASSERT((boost::is_same< typename geometry_traits<T>::tag, typename tag::point3D>));
                if(!inheriter_base::drag_goal) {
                    inheriter_base::drag_goal = this->createGeometry3D(point);
                    inheriter_base::drag_constraint = this->template createConstraint3D<Distance3D>(inheriter_base::drag_point, inheriter_base::drag_goal, 0);
                }
                inheriter_base::drag_goal->set(point);
                this->solve();
            };
            void finishPointDrag() {
                //TODO:remove constraints and drag goal
                inheriter_base::drag_goal.reset();
                inheriter_base::drag_constraint.reset();
            };
        };

        struct inheriter_id : public inheriter_base {

        protected:
            using inheriter_base::m_this;
            using inheriter_base::process_constraint;

        public:
            template<typename T>
            Geom createGeometry3D(T geom, Identifier id) {
                Geom g(new Geometry3D(geom, * ((Sys*) this)));
                fusion::vector<LocalVertex, GlobalVertex> res = m_this->m_cluster.addVertex();
                m_this->m_cluster.template setObject<Geometry3D> (fusion::at_c<0> (res), g);
                g->template setProperty<vertex_prop>(fusion::at_c<1>(res));
                m_this->template objectVector<Geometry3D>().push_back(g);
                g->setIdentifier(id);
                return g;
            };

            template< template<typename,typename,typename> class T >
            Cons createConstraint3D(Identifier id, Geom first, Geom second) {
                Cons c(new Constraint3D(* ((Sys*) this), first, second));
                c->template set<T>(id);
                process_constraint(c, first, second);
                return c;
            };
            template< template<typename,typename,typename> class T, typename Value1 >
            Cons createConstraint3D(Identifier id, Geom first, Geom second, Value1 v1) {
                Cons c(new Constraint3D(* ((Sys*) this), first, second));
                c->template set<T>(id, v1);
                process_constraint(c, first, second);
                return c;
            };
            template< template<typename,typename,typename> class T, typename Value1, typename Value2 >
            Cons createConstraint3D(Identifier id, Geom first, Geom second, Value1 v1, Value2 v2) {
                Cons c(new Constraint3D(* ((Sys*) this), first, second));
                c->template set<T>(id, v1,v2);
                process_constraint(c, first, second);
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

                inheriter_base::drag_point = getGeometry3D(id);
                inheriter_base::drag_goal.reset();
            };

            template<typename T>
            void pointDrag(T point) {
                BOOST_MPL_ASSERT((boost::is_same< typename geometry_traits<T>::tag, typename tag::point3D>));
                if(!inheriter_base::drag_goal) {
                    inheriter_base::drag_goal = this->createGeometry3D(point, "drag_goal");
                    inheriter_base::drag_constraint = this->template createConstraint3D<Fix3D>("drag_constraint", inheriter_base::drag_point, inheriter_base::drag_goal, 0);
                }
                inheriter_base::drag_goal->set(point, "drag_goal");
                //inheriter_base::drag_goal->m_parameterCount=0;
                ((Sys*) this)->solve();
            };
            void finishPointDrag() {
                //TODO:remove constraints and drag goal
                inheriter_base::drag_goal.reset();
                inheriter_base::drag_constraint.reset();
            };
        };

        struct inheriter : public mpl::if_<boost::is_same<Identifier, No_Identifier>, inheriter_noid, inheriter_id>::type {};



        struct math_prop {
            typedef cluster_property kind;
            typedef details::ClusterMath<Sys> type;
        };
        struct gmap_prop {
            typedef cluster_property kind;
            typedef std::vector<Geom> type;
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

        typedef mpl::vector<vertex_prop, edge_prop, math_prop, gmap_prop, fix_prop>  properties;

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

}

#endif //GCM_GEOMETRY3D_H




