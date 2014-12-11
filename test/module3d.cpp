/*
    openDCM, dimensional constraint manager
    Copyright (C) 2012  Stefan Troeger <stefantroeger@gmx.net>

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

#include <boost/test/unit_test.hpp>
#include <boost/exception/get_error_info.hpp>

#include <Eigen/Core>
#include "opendcm/core.hpp"
#include "opendcm/module3d.hpp"
#include "derivativetest.hpp"


typedef dcm::Eigen3Kernel<double> K;

template<typename Kernel, bool Map = true>
struct TDirection3 : public dcm::geometry::Geometry<Kernel, Map,
        dcm::geometry::storage::Vector<3>> {

    typedef typename Kernel::Scalar Scalar;
    using dcm::geometry::Geometry<Kernel, Map, dcm::geometry::storage::Vector<3>>::m_storage;

    auto direction() -> decltype(fusion::at_c<0>(m_storage)) {
        return fusion::at_c<0>(m_storage);
    };

    TDirection3<Kernel, Map>& transform(const Eigen::Transform<Scalar, 3, Eigen::AffineCompact>& t) {
        direction() = t.rotation()*direction();
        return *this;
    };

    TDirection3<Kernel, Map>  transformed(const Eigen::Transform<Scalar, 3, Eigen::AffineCompact>& t) {
        TDirection3<Kernel, Map> copy(*this);
        copy.transform(t);
        return copy;
    };
};

typedef Eigen::Matrix<double, 3,1> Vector3;


typedef dcm::geometry::adaptor<TDirection3> Direction3;
template<>
struct dcm::geometry_traits<Vector3> {
    typedef Direction3 type;
    typedef dcm::modell::CartesianDirection modell;
    typedef dcm::accessor::OrderdBracket    accessor;
};

typedef dcm::System<dcm::Module3D<Vector3>> System;

BOOST_AUTO_TEST_SUITE(Module3D_test_suit);

BOOST_AUTO_TEST_CASE(cluster) {

    typedef dcm::numeric::Cluster3<K>::ParameterIterator  cParIt;
    typedef dcm::numeric::Cluster3<K>::Derivative         cDer;
    typedef dcm::numeric::Cluster3Geometry<K, TDirection3>::Derivative clDer;

    try {

        dcm::numeric::LinearSystem<K> sys(10,10);
        dcm::numeric::Cluster3<K> cluster;

        cluster.init(sys);
        cluster.recalculate();

        BOOST_CHECK(cluster.parameters().size()==6);
        BOOST_CHECK(cluster.derivatives().size()==6);

        dcm::numeric::Cluster3Geometry<K, TDirection3> clGeom;
        clGeom.init(sys);

        BOOST_CHECK(clGeom.parameters().size()==0);
        BOOST_CHECK(clGeom.derivatives().size()==0);

        clGeom.setBaseGeometry(&cluster);
        clGeom.recalculate();

        BOOST_CHECK(clGeom.parameters().size()==6);
        BOOST_CHECK(clGeom.derivatives().size()==6);

        cluster.addClusterGeometry(&clGeom);
        cluster.recalculate();

        for(clDer& der : clGeom.derivatives())
            BOOST_CHECK(der.second.Value != nullptr);

        //let's test the derivatives and see if we calculate them correct for 0 values
        DerivativeTest::isCorrect(clGeom, 
               [&]() {
                   cluster.recalculate();
                   clGeom.recalculate();
               }
        );
        
        //and check the derivatives for an arbitrary value
        clGeom.direction() << 1,2,3;
        clGeom.direction().normalize();
        DerivativeTest::isCorrect(clGeom, 
               [&]() {
                   cluster.recalculate();
                   clGeom.recalculate();
               }
        );


    }
    catch(boost::exception& x) {
        BOOST_FAIL(*boost::get_error_info<dcm::error_message>(x));
    }
    catch(std::exception& x) {
        BOOST_FAIL("Unknown exception");
    }
};

BOOST_AUTO_TEST_CASE(geometry) {
    
    Vector3 v;
    v<<1,2,3;
    
    System s;
    std::shared_ptr<System::Geometry3D> g = s.addGeometry3D(v);
    
    BOOST_CHECK(g->holdsType());
    BOOST_CHECK(g->holdsGeometry());
    BOOST_CHECK(g->holdsGeometryType<Vector3>());
    BOOST_CHECK(g->holdsGeometryType<Direction3>());
    
    //Vector3 v2 = dcm::get<Vector3>(g);
    //BOOST_CHECK(v2.isApprox(v));
    
    Vector3 v3 = g->get<Vector3>();
    BOOST_CHECK(v3.isApprox(v));
};

BOOST_AUTO_TEST_SUITE_END();
