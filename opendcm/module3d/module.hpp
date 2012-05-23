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

#ifndef DCM_MODULE_3D_H
#define DCM_MODULE_3D_H

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

#include "opendcm/core/object.hpp"
#include "opendcm/core/clustergraph.hpp"
#include "opendcm/core/sheduler.hpp"
#include "opendcm/core/traits.hpp"
#include "opendcm/core/geometry.hpp"
#include "geometry.hpp"
#include "constraint.hpp"
#include "dof.hpp"


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
    typename Kernel::Vector3Map m_translation;
    typename Kernel::Quaternion	m_quaternion;
    typename Kernel::Vector3	m_original_translation;
    typename Kernel::Vector3Map	m_normQ;

public:
    ClusterMath() : m_normQ(NULL), m_translation(NULL) {};

    void setRotationMap(typename Kernel::Matrix3Map& map, typename Kernel::Matrix39Map& diffmap) {
        new(&map) typename Kernel::Matrix3Map(&m_rotation(0,0),3,3);
        new(&diffmap) typename Kernel::Matrix39Map(&m_diffrot(0,0));
    };
    void setTranslationMap(typename Kernel::Vector3Map& map) {
        new(&map) typename Kernel::Vector3Map(&m_translation(0));
    };
    typename Kernel::Vector3Map& getNormQuaternionMap() {
        return m_normQ;
    };
    typename Kernel::Vector3Map& getTranslationMap() {
        return m_translation;
    };
    void initMaps() {

        //nQ = a,b,c ; n = ||nQ»» ;  Q = x,y,z,w = a/n,b/n,c/n,n
        // --> w = n; a = x*n = x*w ...
        m_normQ(0) = m_quaternion.x()*m_quaternion.w();
        m_normQ(1) = m_quaternion.y()*m_quaternion.w();
        m_normQ(2) = m_quaternion.z()*m_quaternion.w();
        m_translation = m_original_translation;
    };

    typename Kernel::Quaternion& getQuaternion() {
        return m_quaternion;
    };
    typename Kernel::Vector3& getTranslation() {
        return m_original_translation;
    };

    void setValuesFromMaps() {
        Scalar norm = m_normQ.norm();
        m_quaternion = typename Kernel::Quaternion(m_normQ(0)/norm, m_normQ(1)/norm, m_normQ(2)/norm, norm);
        m_quaternion.normalize();
        m_original_translation = m_translation;
    };

    void recalculate() {

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
        typedef boost::shared_ptr<Constraint3D> Cons;

        typedef mpl::map< mpl::pair<reset, boost::function<void (Geom) > > >  GeomSignal;

        struct MES  : public system_traits<Sys>::Kernel::MappedEquationSystem {

            typedef typename system_traits<Sys>::Cluster Cluster;
            Cluster& m_cluster;

            MES(Cluster& cl, int par, int eqn) : system_traits<Sys>::Kernel::MappedEquationSystem(par, eqn),
                m_cluster(cl) {};

            virtual void recalculate() {

                //first calculate all clusters
                typedef typename Cluster::cluster_iterator citer;
                std::pair<citer, citer> cit = m_cluster.clusters();
                for(; cit.first != cit.second; cit.first++) {

                    (*cit.first).second->template getClusterProperty<math_prop>().recalculate();

                    //now with the new rotation matrix we calculate all geometries in that cluster
                    std::vector<Geom>& vec = (*cit.first).second->template getClusterProperty<gmap_prop>();
                    typedef typename std::vector<Geom>::iterator iter;

                    for(iter it = vec.begin(); it != vec.end(); it++)
                        (*it)->recalculate();

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
                            (*oit.first)->calculate();
                    }
                }
            }
        };

        struct SystemSolver : public Job<Sys> {

            typedef typename system_traits<Sys>::Cluster Cluster;
            typedef typename system_traits<Sys>::Kernel Kernel;

            SystemSolver() {
                Job<Sys>::priority = 1000;
            };

            virtual void execute(Sys& sys) {

                std::cout<<"Execute system solver"<<std::endl;
                solveCluster(sys.m_cluster);
            };

            void solveCluster(Cluster& cluster) {

                uint parameters=0, constraints=0, offset=0, equation=0;

                //get the ammount of parameters and constraint equations we need
                typedef typename boost::graph_traits<Cluster>::vertex_iterator iter;
                std::pair<iter, iter>  it = boost::vertices(cluster);
                for(; it.first != it.second; it.first++) {

                    if(cluster.isCluster(*it.first)) parameters += 6;
                    else {
                        parameters += cluster.template getObject<Geometry3D>(*it.first)->m_parameterCount;
                    };
                }

                typedef typename boost::graph_traits<Cluster>::edge_iterator e_iter;
                std::pair<e_iter, e_iter>  e_it = boost::edges(cluster);
                for(; e_it.first != e_it.second; e_it.first++)
                    constraints += cluster.getGlobalEdgeCount(*e_it.first);

                //initialise the system with now known size
                MES mes(cluster, parameters, constraints);

                //iterate all geometrys again and set the needed maps
                it = boost::vertices(cluster);
                for(; it.first != it.second; it.first++) {

                    if(cluster.isCluster(*it.first)) {
                        //set norm Quaternion as map to the parameter vector
                        Cluster& c = cluster.getVertexCluster(*it.first);
                        details::ClusterMath<Sys>& cm =  c.template getClusterProperty<math_prop>();
                        mes.setParameterMap(offset, cm.getNormQuaternionMap());
                        //set translation as map to the parameter vector
                        mes.setParameterMap(offset+3, cm.getTranslationMap());
                        //write initail values to the parameter maps
                        cm.initMaps();

                        //map all geometrie within that cluster to it's rotation matrix
                        std::vector<Geom>& vec = cluster.template getClusterProperty<gmap_prop>();
                        vec.clear();
                        mapClusterDownstreamGeometry(c, cm, vec);

                        offset += 6;
                    } else {
                        Geom g = cluster.template getObject<Geometry3D>(*it.first);
                        mes.setParameterMap(offset, g->m_parameterCount, g->getParameterMap());
                        g->m_parameterOffset = offset;
                        //init the parametermap with initial values
                        g->initMap();
                        offset += g->m_parameterCount;
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
                            mes.setResidualMap(equation, c->m_residual);
                            mes.setJacobiMap(equation, c->first->m_parameterOffset, c->first->m_parameterCount, c->m_diffFirst);
                            mes.setJacobiMap(equation, c->second->m_parameterOffset, c->second->m_parameterCount, c->m_diffSecond);
                            equation++;
                        }
                        //TODO: else throw (as every global edge was counted as one equation)
                    }
                }

                std::cout<<"Residual bevore solving: "<<mes.Residual.norm()<<std::endl;

                //now it's time to solve
                Kernel::solve(mes);

                std::cout<<"Residual after solving: "<<mes.Residual.norm()<<std::endl;
            };

            void mapClusterDownstreamGeometry(Cluster& cluster, details::ClusterMath<Sys>& cm, std::vector<Geom>& vec) {
                //all geometry within that cluster needs to be mapped to the provided rotation matrix (in cm)

                //get all vertices and map the geometries if existend
                typedef typename boost::graph_traits<Cluster>::vertex_iterator iter;
                std::pair<iter, iter>  it = boost::vertices(cluster);
                for(; it.first != it.second; it.first++) {
                    Geom g = cluster.template getObject<Geometry3D>(*it.first);
                    if(g) {
                        //allow iteration over all maped geometries
                        vec.push_back(g);
                        //map rotation and diffrotation from cluster to geometry
                        cm.setRotationMap(g->getRotationMap(), g->getDiffRotationMap());
                        //map translation from cluster to geometry
                        cm.setTranslationMap(g->getTranslationMap());

                        g->setClusterMode(true);
                    }
                }

                //go downstream and map
                typedef typename Cluster::cluster_iterator citer;
                std::pair<citer, citer> cit = cluster.clusters();
                for(; cit.first != cit.second; cit.first++)
                    mapClusterDownstreamGeometry(*(*cit.first).second, cm, vec);

            };

        };

        class Geometry3D : public Object<Sys, Geometry3D, GeomSignal > {
            typedef typename boost::make_variant_over< Typelist >::type Variant;
            typedef Object<Sys, Geometry3D, GeomSignal> base;
            typedef typename system_traits<Sys>::Kernel Kernel;
            typedef typename Kernel::number_type Scalar;
            typedef typename Kernel::DynStride DS;

        public:
            template<typename T>
            Geometry3D(T geometry, Sys& system) : base(system), m_isInCluster(false),
                m_geometry(geometry), m_rotation(NULL), m_parameter(NULL,0,DS(0,0)),
                m_diffrot(NULL), m_translation(NULL)  {

                init<T>(geometry);
            };


            template<typename T>
            void set(T geometry) {
                m_geometry = geometry;
                init<T>(geometry);
                base::template emitSignal<reset> (base::shared_from_this());
            };

        protected:
            Variant m_geometry;
            int     m_parameterCount, m_parameterOffset, m_transformations;
            bool    m_isInCluster;
            typename Sys::Kernel::Vector      m_original, m_value; //original and rotated parameters
            typename Sys::Kernel::Matrix      m_diffparam; //gradient vectors combined as matrix when in cluster
            typename Sys::Kernel::VectorMap   m_parameter; //map to the parameters in the solver
            typename Sys::Kernel::Vector3Map  m_translation; //map to the cluster translation
            typename Sys::Kernel::Matrix3Map  m_rotation; //map to the cluster rotation
            typename Sys::Kernel::Matrix39Map m_diffrot; //map to the gradient rotations

            template<typename T>
            void init(T& t) {
                m_parameterCount = geometry_traits<T>::tag::parameters::value;
                m_transformations = geometry_traits<T>::tag::transformations::value;

                m_original.resize(m_parameterCount);
                m_value.resize(m_parameterCount);

                m_diffparam.resize(m_parameterCount,6);
                //the translation gradient is always constant
                for(int i=0; i!=m_transformations; i++)
                    m_diffparam.block(i*3,3,3,3).setIdentity();
                //the not transformed parameters (constant when in cluster) have gradient 0
                int ntp = m_parameterCount-m_transformations*3; //not transformed parameter
                m_diffparam.block(m_parameterCount-ntp, 0, ntp, 6).setZero();

                (typename geometry_traits<T>::modell()).template extract<Scalar,
                typename geometry_traits<T>::accessor >(t, m_original);
            }

            typename Sys::Kernel::VectorMap& getParameterMap() {
                m_isInCluster = false;
                //TODO: write initial values into m_parameter
                return m_parameter;
            }
            typename Sys::Kernel::Matrix3Map&  getRotationMap() {
                return m_rotation;
            };
            typename Sys::Kernel::Matrix39Map& getDiffRotationMap() {
                return m_diffrot;
            };
            typename Sys::Kernel::Vector3Map&  getTranslationMap() {
                return m_translation;
            };
            void initMap() {
                m_parameter = m_original;
            };

            void setClusterMode(bool iscluster) {
                m_isInCluster = iscluster;
                if(iscluster) {
                    //we are in cluster, therfore the parameter map should not point to a solver value but to
                    //the rotated original value;
                    new(&m_parameter) typename Sys::Kernel::VectorMap(&m_value(0), m_parameterCount);
                };
            }
            bool getClusterMode() {
                return m_isInCluster;
            };

            void recalculate() {
                if(!m_isInCluster) return;

                for(int i=0; i!=m_transformations; i++) {
                    //first transform the original to the transformed value
                    m_value.block(i*3,0,3,0) = m_rotation*m_original.block(i*3,0,3,0);

                    //now calculate the gradient vectors and add them to diffparam
                    m_diffparam.block(i*3,0,3,1) = m_diffrot.block(0,0,3,3) * m_original(i*3,3);
                    m_diffparam.block(i*3,1,3,1) = m_diffrot.block(0,3,3,3) * m_original(i*3,3);
                    m_diffparam.block(i*3,2,3,1) = m_diffrot.block(0,6,3,3) * m_original(i*3,3);
                }
            }

            friend class SystemSolver;
            friend class Constraint3D;
            friend struct MES;
        };

        //type erasure container for constraints
        class Constraint3D : public Object<Sys, Constraint3D, mpl::map<> > {

            typedef typename system_traits<Sys>::Kernel Kernel;
            typedef typename Kernel::number_type Scalar;
            typedef typename Kernel::DynStride DS;

        public:
            Constraint3D(Sys& system, Geom f, Geom s) : Object<Sys, Constraint3D, mpl::map<> > (system),
                first(f), second(s), content(0), m_diffFirst(NULL,0,DS(0,0)), m_diffSecond(NULL,0,DS(0,0)),
                m_residual(NULL,0,DS(0,0))	{

                cf = first->template connectSignal<reset> (boost::bind(&Constraint3D::geometryReset, this, _1));
                cs = second->template connectSignal<reset> (boost::bind(&Constraint3D::geometryReset, this, _1));
            };

            ~Constraint3D()  {
                delete content;
                first->template disconnectSignal<reset>(cf);
                second->template disconnectSignal<reset>(cs);
            }

            template< template<typename,typename,typename> class T>
            void setType() {
                creator<T> creator;
                boost::apply_visitor(creator, first->m_geometry, second->m_geometry);
                content = creator.p;
		if(creator.need_swap) first.swap(second);
            };

        protected:
            template<typename T> void pretty(T t) {
                std::cout<<__PRETTY_FUNCTION__<<std::endl;
            }

            Scalar calculate() {

                //first the residual (operator= doeas not work with scalars)
                m_residual << content->calculate(first->m_parameter, second->m_parameter);

                //now see which way we should calculate the gradient (may be diffrent for both geometries)
                if(first->getClusterMode()) {
                    //cluster mode, so we do a full calculation with all 6 diffparam vectors
                    std::cout<<"test"<<std::endl;
                    pretty(m_diffFirst.block(0,0,1,1));
                } else {
                    //not in cluster, so allow the constraint to optimize the gradient calculation
                    content->calculateFirstGradient(first->m_parameter, second->m_parameter, m_diffFirst);
                }
                if(second->getClusterMode()) {
                    //cluster mode, so we do a full calculation with all 6 diffparam vectors
                    std::cout<<"test"<<std::endl;
                    pretty(m_diffFirst.block(0,0,1,1));
                } else {
                    //not in cluster, so allow the constraint to optimize the gradient calculation
                    content->calculateFirstGradient(first->m_parameter, second->m_parameter, m_diffSecond);
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
                virtual Scalar calculateFirstFullGradient(typename Kernel::VectorMap& param1,
                        typename Kernel::VectorMap& param2,
                        typename Kernel::VectorMap& diffparam) = 0;
                virtual Scalar calculateSecondFullGradient(typename Kernel::VectorMap& param1,
                        typename Kernel::VectorMap& param2,
                        typename Kernel::VectorMap& diffparam) = 0;
                virtual void calculateFirstGradient(typename Kernel::VectorMap& param1,
                                                    typename Kernel::VectorMap& param2,
                                                    typename Kernel::VectorMap& grad) = 0;
                virtual void calculateSecondGradient(typename Kernel::VectorMap& param1,
                                                     typename Kernel::VectorMap& param2,
                                                     typename Kernel::VectorMap& grad) = 0;
            };

            template< template<typename, typename, typename> class T1, typename T2, typename T3>
            struct holder : public placeholder  {

                holder(const T1<Kernel, T2, T3> & value)
                    : held(value)   {}

                virtual Scalar calculate(typename Kernel::VectorMap& f, typename Kernel::VectorMap& s) {
                    return held.calculate(f,s);
                };
                virtual Scalar calculateFirstFullGradient(typename Kernel::VectorMap& param1,
                        typename Kernel::VectorMap& param2,
                        typename Kernel::VectorMap& diffparam) {
                    return held.calculateSecondFullGradient(param1, param2, diffparam);
                };
                virtual Scalar calculateSecondFullGradient(typename Kernel::VectorMap& param1,
                        typename Kernel::VectorMap& param2,
                        typename Kernel::VectorMap& diffparam) {
                    return held.calculateSecondFullGradient(param1, param2, diffparam);
                };
                virtual void calculateFirstGradient(typename Kernel::VectorMap& param1,
                                                    typename Kernel::VectorMap& param2,
                                                    typename Kernel::VectorMap& grad) {
                    held.calculateFirstGradient(param1, param2, grad);
                };
                virtual void calculateSecondGradient(typename Kernel::VectorMap& param1,
                                                     typename Kernel::VectorMap& param2,
                                                     typename Kernel::VectorMap& grad) {
                    held.calculateSecondGradient(param1, param2, grad);
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

            placeholder* content;
            Geom first, second;
            Connection cf, cs;

            typename Kernel::VectorMap m_diffFirst, m_diffSecond, m_residual;

            friend class SystemSolver;
            friend class MES;

        };

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
                m_this->template objectVector<Geometry3D>().push_back(g);
                return g;
            };

            template< template<typename,typename,typename> class T >
            Cons createConstraint3D(Geom first, Geom second) {

                Cons c(new Constraint3D(* ((Sys*) this), first, second));
                c->template setType<T>();
                fusion::vector<LocalEdge, GlobalEdge, bool, bool> res;
                res = m_this->m_cluster.addEdge(first->template getProperty<vertex_prop>(),
                                                second->template getProperty<vertex_prop>());
                m_this->m_cluster.template setObject<Constraint3D> (fusion::at_c<0> (res), c);
                c->template setProperty<edge_prop>(fusion::at_c<1>(res));
                m_this->template objectVector<Constraint3D>().push_back(c);
                return c;
            };

        private:
            Sys* m_this;
        };


        struct math_prop {
            typedef cluster_property kind;
            typedef details::ClusterMath<Sys> type;
        };
        struct gmap_prop {
            typedef cluster_property kind;
            typedef std::vector<Geom> type;
        };
        struct vertex_prop {
            typedef Geometry3D kind;
            typedef GlobalVertex type;
        };
        struct edge_prop {
            typedef Constraint3D kind;
            typedef GlobalEdge type;
        };

        typedef mpl::vector<vertex_prop, edge_prop, math_prop, gmap_prop>  properties;

        static void system_init(Sys& sys) {
            std::cout<<"add solver job"<<std::endl;
            sys.m_sheduler.addProcessJob(new SystemSolver());
        };
    };
};

}

#endif //DCM_GEOMETRY3D_H
