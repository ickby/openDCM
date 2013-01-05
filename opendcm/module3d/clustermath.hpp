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

#ifndef GCM_CLUSTERMATH_H
#define GCM_CLUSTERMATH_H

#include <vector>
#include <cmath>

#define maxfak 1.2
#define minfak 0.8

namespace dcm {
namespace details {

enum Scalemode {
    one,
    two,
    three,
    multiple_inrange,
    multiple_outrange
};

template<typename Sys>
struct ClusterMath {

public:
    typedef typename system_traits<Sys>::Kernel Kernel;
    typedef typename system_traits<Sys>::Cluster Cluster;
    typedef typename system_traits<Sys>::template getModule<m3d>::type module3d;
    typedef typename module3d::Geometry3D Geometry3D;
    typedef boost::shared_ptr<Geometry3D> Geom;
    typedef typename module3d::type::math_prop math_prop;
    typedef typename module3d::type::fix_prop fix_prop;

    typedef typename Kernel::number_type Scalar;

    typename Kernel::Matrix3	 m_rotation;
    typename Kernel::Matrix39 	 m_diffrot;
    typename Kernel::Transform3D m_transform;
    typename Kernel::Vector3Map	 m_normQ;
    typename Kernel::Quaternion  m_resetQuaternion;

    int m_offset;
    bool init, reset;
    std::vector<Geom> m_geometry;

    typename Kernel::Vector3Map m_translation;
    //shift scale stuff
    typename Kernel::Vector3 midpoint, m_shift, scale_dir, maxm, minm, max, fixtrans;
    Scalar m_scale;
    Scalemode mode;

public:
    ClusterMath() : m_normQ(NULL), m_translation(NULL), init(false),
        m_resetQuaternion(1,1,1,1), reset(false) {

        m_resetQuaternion.normalize();
        m_shift.setZero();
        m_scale = 1.;
    };

    void setParameterOffset(int offset) {
        m_offset = offset;
    };
    int getParameterOffset() {
        return m_offset;
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

        //if we start with unit quaternion we need to reset the geometry to something else
        const typename Kernel::Quaternion& m_quaternion = m_transform.rotation();
        if(m_quaternion.w() < 1.) {
            Scalar s = std::acos(m_quaternion.w())/std::sin(std::acos(m_quaternion.w()));
            m_normQ = m_quaternion.vec()*s;
        } else {
            m_normQ.setZero();
        }

        m_translation = m_transform.translation().vector();
        init = true;
        m_scale = 1.;
        midpoint.setZero();
        m_shift.setZero();
    };
    void initFixMaps() {
        //when fixed no maps exist
        m_diffrot.setZero();
        m_rotation = m_transform.rotation().toRotationMatrix();
        new(&m_translation) typename Kernel::Vector3Map(&fixtrans(0));
	m_translation = m_transform.translation().vector();
        init = true;
        m_scale = 1.;
        midpoint.setZero();
        m_shift.setZero();
    };

    typename Kernel::Transform3D& getTransform() {
        return m_transform;
    };
    void setShift(typename Kernel::Vector3 s) {
        m_shift = s;
        //we remove shift from the local geometries, therefore we have to add it here
        //to not change the global position
        if(init) 
             m_translation += m_transform*m_shift;

    };
    void setScale(Scalar s) {
        m_scale = s;
        if(init)m_translation *= s;
    };

    void finishCalculation() {
        const Scalar norm = m_normQ.norm();
        const Scalar fac = std::sin(norm)/norm;
        typename Kernel::Quaternion m_quaternion = typename Kernel::Quaternion(std::cos(norm), m_normQ(0)*fac, m_normQ(1)*fac, m_normQ(2)*fac);
        m_quaternion.normalize();

        //if the rotation was resetet we have to calc the real quaternion now
        if(reset)
            resetClusterRotation(m_quaternion);


        typename Kernel::Vector3 m_original_translation = m_translation/m_scale - m_quaternion._transformVector(m_shift);

        //needed to allow a correct global calculation in cluster geometries after this finish
        m_shift.setZero();
        m_rotation = m_quaternion.toRotationMatrix();
        m_translation = m_original_translation;

        init=false;
	
	m_transform = m_quaternion;
	m_transform *= typename Kernel::Transform3D::Translation(m_original_translation);
	m_transform *= typename Kernel::Transform3D::Scaling(m_scale);
    };
    void finishFixCalculation() {
        //needed to allow a correct global calculation in cluster geometries after this finish
        m_shift.setZero();
        m_translation /= m_scale;
    }

    void resetClusterRotation(typename Kernel::Quaternion& Q) {

        typename Kernel::Matrix3 rot;
        const typename Kernel::Vector3 trans(0,0,0);
        if(!reset) {
            Q = Q*m_resetQuaternion.conjugate();
            rot = m_resetQuaternion.toRotationMatrix();
            reset = true;
        } else {
            Q = Q*m_resetQuaternion;
            rot = m_resetQuaternion.conjugate().toRotationMatrix();
            reset = false;
        }

        //apply the needed transformation to all geometries local values
        typedef typename std::vector<Geom>::iterator iter;
        for(iter it = m_geometry.begin(); it != m_geometry.end(); it++) {
            (*it)->transformLocal(rot, trans);
        };
        //also transform the shift (as it is calculated from the local values which are transformed now)
        m_shift = rot*m_shift;

        if(Q.w() < 1.) {
            //set the normQ map to the new quaternion value
            const Scalar s = std::acos(Q.w())/std::sin(std::acos(Q.w()));
            m_normQ = Q.vec()*s;
        } else m_normQ.setZero();
    };

    void recalculate() {

        typename Kernel::Quaternion Q;
        Scalar norm = m_normQ.norm();

        //prevent the quaternion to be the idendity form
        if(norm == 0) {
            Q.setIdentity();
            resetClusterRotation(Q);
            norm = m_normQ.norm();
        }

        //get the Quaternion for the norm quaternion form and calculate the rotation matrix
        const Scalar fac = std::sin(norm)/norm;
        Q = typename Kernel::Quaternion(std::cos(norm), m_normQ(0)*fac, m_normQ(1)*fac, m_normQ(2)*fac);
        Q.normalize(); //not needed, just to avoid rounding errors

        //also being too close to the identity sucks
        if(norm < 0.1)
            resetClusterRotation(Q);

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

    void addGeometry(Geom g) {
        m_geometry.push_back(g);
    };
    void clearGeometry() {
        m_geometry.clear();
    };
    std::vector<Geom>& getGeometry() {
        return m_geometry;
    };


    void mapClusterDownstreamGeometry(Cluster& cluster,
                                      typename Kernel::Transform3D& t,
                                      details::ClusterMath<Sys>& cm
                                     ) {
        //all geometry within that cluster needs to be mapped to the provided rotation matrix (in cm)
        //also the geometries toplocal value needs to be set so that it matches this cm
        typename Kernel::Transform3D nt = t*cluster.template getClusterProperty<math_prop>().getTransform();
	typename Kernel::Transform3D trans = nt.inverse();
		
        //get all vertices and map the geometries if existend
        typedef typename boost::graph_traits<Cluster>::vertex_iterator iter;
        std::pair<iter, iter>  it = boost::vertices(cluster);
        for(; it.first != it.second; it.first++) {
            Geom g = cluster.template getObject<Geometry3D>(*it.first);
            if(g) {
                //allow iteration over all maped geometries
                cm.addGeometry(g);
                //map rotation and diffrotation from cluster to geometry
                cm.setRotationMap(g->getRotationMap(), g->getDiffRotationMap());
                //map translation from cluster to geometry
                cm.setTranslationMap(g->getTranslationMap());
                //map shift from cluster to geometry
                cm.setShiftMap(g->getShiftMap());
                //set the offsets so that geometry knows where it is in the parameter map
                g->m_offset = cm.getParameterOffset();
                //calculate the appropriate local values
                g->transform(trans);

                //position and offset of the parameters must be set to the clusters values
                g->setClusterMode(true, cluster.template getClusterProperty<fix_prop>());
            }
        }

        //go downstream and map
        typedef typename Cluster::cluster_iterator citer;
        std::pair<citer, citer> cit = cluster.clusters();
        for(; cit.first != cit.second; cit.first++)
            mapClusterDownstreamGeometry(*(*cit.first).second, nt, cm);
        //TODO: if one subcluster is fixed the hole cluster should be too, as there are no
        //	dof's remaining between parts and so nothing can be moved when one part is fixed.

    };

    /*Calculate the scale of the cluster. Therefore the midpoint is calculated and the scale is
     * defined as the max distance between the midpoint and the points.
    */
    Scalar calculateClusterScale() {

        if(m_geometry.empty()) assert(false); //TODO: Throw
        else if(m_geometry.size() == 1) {
            return calcOnePoint(m_geometry[0]->getPoint());
        } else if(m_geometry.size() == 2) {
            const typename Kernel::Vector3 p1 = m_geometry[0]->getPoint();
            const typename Kernel::Vector3 p2 = m_geometry[1]->getPoint();

            if(Kernel::isSame((p1-p2).norm(), 0.))
                return calcOnePoint(p1);

            return calcTwoPoints(p1, p2);
        } else if(m_geometry.size() == 3) {

            const typename Kernel::Vector3 p1 = m_geometry[0]->getPoint();
            const typename Kernel::Vector3 p2 = m_geometry[1]->getPoint();
            const typename Kernel::Vector3 p3 = m_geometry[2]->getPoint();

            const typename Kernel::Vector3 d = p2-p1;
            const typename Kernel::Vector3 e = p3-p1;

            if(Kernel::isSame(d.norm(), 0.)) {

                if(Kernel::isSame(e.norm(), 0.))
                    return calcOnePoint(p1);

                return calcTwoPoints(p1, p3);
            } else if(Kernel::isSame(e.norm(), 0.)) {
                return calcTwoPoints(p1, p2);
            } else if(!Kernel::isSame((d/d.norm() - e/e.norm()).norm(), 0.) &&
                      !Kernel::isSame((d/d.norm() + e/e.norm()).norm(), 0.)) {
                return calcThreePoints(p1, p2, p3);
            }
            //three points on a line need to be treaded as multiple points
        }

        //more than 3 points dont have a exakt solution. we search for a midpoint from which all points
        //are at least maxfak*scale away, but not closer than minfak*scale

        //get the bonding box to get the center of points
        Scalar xmin=1e10, xmax=1e-10, ymin=1e10, ymax=1e-10, zmin=1e10, zmax=1e-10;
        typedef typename std::vector<Geom>::iterator iter;
        for(iter it = m_geometry.begin(); it != m_geometry.end(); it++) {
            typename Kernel::Vector3 v = (*it)->getPoint();
            xmin = (v(0)<xmin) ? v(0) : xmin;
            xmax = (v(0)<xmax) ? xmax : v(0);
            ymin = (v(1)<ymin) ? v(1) : ymin;
            ymax = (v(1)<ymax) ? ymax : v(1);
            zmin = (v(2)<zmin) ? v(2) : zmin;
            zmax = (v(2)<zmax) ? zmax : v(2);
        };
        //now calculate the midpoint
        midpoint << xmin+xmax, ymin+ymax, zmin+zmax;
        midpoint /= 2.;


        //get the scale direction an the resulting nearest point indexes
        double xh = xmax-xmin;
        double yh = ymax-ymin;
        double zh = zmax-zmin;
        int i1, i2, i3;
        if((xh<=yh) && (xh<=zh)) {
            i1=1;
            i2=2;
            i3=0;
        } else if((yh<xh) && (yh<zh)) {
            i1=0;
            i2=2;
            i3=1;
        } else {
            i1=0;
            i2=1;
            i3=2;
        }
        scale_dir.setZero();
        scale_dir(i3) = 1;
        max = Eigen::Vector3d(xmin,ymin,zmin);
        m_scale = (midpoint-max).norm()/maxfak;
        mode = multiple_inrange;

        maxm = max-midpoint;
        Scalar minscale = 1e10;

        //get the closest point
        for(iter it = m_geometry.begin(); it != m_geometry.end(); it++) {

            const Eigen::Vector3d point = (*it)->getPoint()-midpoint;
            if(point.norm()<minfak*m_scale) {

                const double h = std::abs(point(i3)-maxm(i3));
                const double k = std::pow(minfak/maxfak,2);
                double q = std::pow(point(i1),2) + std::pow(point(i2),2);
                q -= (std::pow(maxm(i1),2) + std::pow(maxm(i2),2) + std::pow(h,2))*k;
                q /= 1.-k;

                const double p = h*k/(1.-k);

                if(std::pow(p,2)<q) assert(false);

                midpoint(i3) += p + std::sqrt(std::pow(p,2)-q);
                maxm = max-midpoint;
                m_scale = maxm.norm()/maxfak;

                mode = multiple_outrange;
                minm = (*it)->getPoint()-midpoint;

                it = m_geometry.begin();
            } else if(point.norm()<minscale) {
                minscale = point.norm();
            }
        }

        if(mode==multiple_inrange) {
            //we are in the range, let's get the perfect balanced scale value
            m_scale = (minscale+maxm.norm())/2.;
        }
        return m_scale;
    };

    void applyClusterScale(Scalar scale, bool isFixed) {

        //when fixed, the geometries never get recalculated. therefore we have to do a calculate now
        //to alow the adoption of the scale. and no shift should been set.
        if(isFixed) {
            //scale needs to be set to have the correct translation value
            setScale(1./scale);
            //now calculate the scaled geometrys
            typedef typename std::vector<Geom>::iterator iter;
            for(iter it = m_geometry.begin(); it != m_geometry.end(); it++) {
                (*it)->recalculate(1./scale);
            };
            return;
        }

        //if this is our scale then just applie the midpoint as shift
        if(Kernel::isSame(scale, m_scale)) {
            setShift(midpoint);
            setScale(1./scale);
            return;
        }

        //if only one point exists we extend the origin-point-line to match the scale
        if(mode==details::one) {
            if(Kernel::isSame(midpoint.norm(),0))
                m_shift << scale, 0, 0;
            else m_shift = midpoint + scale*scale_dir;

            setShift(m_shift);
            setScale(1./scale);
        }
        //two and three points form a rectangular triangle, so same procedure
        else if(mode==details::two || mode==details::three) {

            setShift(midpoint+scale_dir*std::sqrt(std::pow(scale,2) - std::pow(m_scale,2)));
            setScale(1./scale);
        }
        //multiple points
        else if(mode==details::multiple_outrange) {

            if(scale_dir(0)) {
                Scalar d = std::pow(maxm(1),2) + std::pow(maxm(2),2);
                Scalar h = std::sqrt(std::pow(maxfak*scale,2)-d);
                midpoint(0) += maxm(0) + h;
            } else if(scale_dir(1)) {
                Scalar d = std::pow(maxm(0),2) + std::pow(maxm(2),2);
                Scalar h = std::sqrt(std::pow(maxfak*scale,2)-d);
                midpoint(1) += maxm(1) + h;
            } else {
                Scalar d = std::pow(maxm(0),2) + std::pow(maxm(1),2);
                Scalar h = std::sqrt(std::pow(maxfak*scale,2)-d);
                midpoint(2) += maxm(2) + h;
            }
            setShift(midpoint);
            setScale(1./scale);
        } else {

            //TODO: it's possible that for this case we get too far away from the outer points.
            //	    The m_scale for "midpoint outside the bounding box" may be bigger than the
            //      scale to applie, so it results in an error.
            //get the closest point
            typedef typename std::vector<Geom>::iterator iter;
            for(iter it = m_geometry.begin(); it != m_geometry.end(); it++) {

                const Eigen::Vector3d point = (*it)->getPoint()-midpoint;
                if(point.norm()<minfak*scale) {

                    if(scale_dir(0)) {
                        Scalar d = std::pow(point(1),2) + std::pow(point(2),2);
                        Scalar h = std::sqrt(std::pow(minfak*scale,2)-d);
                        midpoint(0) += point(0) + h;
                    } else if(scale_dir(1)) {
                        Scalar d = std::pow(point(0),2) + std::pow(point(2),2);
                        Scalar h = std::sqrt(std::pow(minfak*scale,2)-d);
                        midpoint(1) += point(1) + h;
                    } else {
                        Scalar d = std::pow(point(0),2) + std::pow(point(1),2);
                        Scalar h = std::sqrt(std::pow(minfak*scale,2)-d);
                        midpoint(2) += point(2) + h;
                    }
                }
            }
            setShift(midpoint);
            setScale(1./scale);
        }
    };

private:
    Scalar calcOnePoint(typename Kernel::Vector3 p) {

        //one point can have every scale when moving the midpoint on the origin - point vector
        midpoint  = p;
        scale_dir = -midpoint;
        scale_dir.normalize();
        mode = details::one;
        m_scale = 0.;
        return 0.;
    };

    Scalar calcTwoPoints(typename Kernel::Vector3 p1, typename Kernel::Vector3 p2) {

        //two points have their minimal scale at the mid position. Scaling perpendicular to this
        //line allows arbitrary scale values. Best is to have the scale dir move towards the origin
        //as good as possible.
        midpoint  = p1+(p2-p1)/2.;
        scale_dir = (p2-p1).cross(midpoint);
        scale_dir = scale_dir.cross(p2-p1);
        if(!Kernel::isSame(scale_dir.norm(),0)) scale_dir.normalize();
        else scale_dir(0) = 1;
        mode = details::two;
        m_scale = (p2-p1).norm()/2.;
        return m_scale;
    };

    Scalar calcThreePoints(typename Kernel::Vector3 p1,
                           typename Kernel::Vector3 p2, typename Kernel::Vector3 p3) {

        //Three points form a triangle with it's minimal scale at the center of it's outer circle.
        //Arbitrary scale values can be achieved by moving perpendicular to the triangle plane.
        typename Kernel::Vector3 d = p2-p1;
        typename Kernel::Vector3 e = p3-p1;

        typename Kernel::Vector3 f = p1+0.5*d;
        typename Kernel::Vector3 g = p1+0.5*e;
        scale_dir = d.cross(e);

        typename Kernel::Matrix3 m;
        m.row(0) = d.transpose();
        m.row(1) = e.transpose();
        m.row(2) = scale_dir.transpose();

        typename Kernel::Vector3 res(d.transpose()*f, e.transpose()*g, scale_dir.transpose()*p1);

        midpoint =  m.colPivHouseholderQr().solve(res);
        scale_dir.normalize();
        mode = details::three;
        m_scale = (midpoint-p1).norm();

        return m_scale;

    }

};

}//details
}//dcm


#endif //GCM_CLUSTERMATH_H




