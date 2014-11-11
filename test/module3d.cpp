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

#include "opendcm/core/defines.hpp"
#include "opendcm/module3d.hpp"
#include "derivativetest.hpp"


using namespace dcm;

typedef Eigen3Kernel<double> K;

template<typename Kernel, bool Map = true>
struct TDirection3 : public geometry::Geometry<Kernel, Map,
            geometry::storage::Vector<3>> {

    typedef typename Kernel::Scalar Scalar;
    using geometry::Geometry<Kernel, Map, geometry::storage::Vector<3>>::m_storage;
    
    auto value() -> decltype(fusion::at_c<0>(m_storage)){
        return fusion::at_c<0>(m_storage);
    };
    
     TDirection3<Kernel, Map>& transform  (const Eigen::Transform<Scalar, 3, Eigen::AffineCompact>& t) {
         value() = t*value();
         return *this;
     };
     
     TDirection3<Kernel, Map>  transformed(const Eigen::Transform<Scalar, 3, Eigen::AffineCompact>& t) {
         TDirection3<Kernel, Map> copy(*this);
         copy.value() = t*value();
         return copy;
     };
};

BOOST_AUTO_TEST_SUITE(Module3D_test_suit);

BOOST_AUTO_TEST_CASE(cluster) {
   
    typedef numeric::Cluster3<K>::ParameterIterator  cParIt;
    typedef numeric::Cluster3<K>::Derivative         cDer;
    
    try {

        numeric::LinearSystem<K> sys(10,10);
        numeric::Cluster3<K> cluster;
        
        cluster.init(sys);      
        cluster.recalculate();
        
        BOOST_CHECK(cluster.parameters().size()==6);
        BOOST_CHECK(cluster.derivatives().size()==6);
                
        int c = 0;
        cParIt it = cluster.parameters().begin();
        for(cDer& der : cluster.derivatives()) {
            if(c>2) {
                BOOST_CHECK(der.first.translation()(c-3)==1);
                BOOST_CHECK(der.first.translation().sum()==1);
                BOOST_CHECK(der.first.rotation().sum()==0);
            } 
            else {
                BOOST_CHECK(der.first.translation().sum()==0);
                BOOST_CHECK(der.first.rotation().sum()!=0);
            }
            BOOST_CHECK(der.second == *(it++));
            BOOST_CHECK(der.second.Value != nullptr);
            ++c;
            
            std::cout<<"translation: "<<der.first.translation().transpose()<<std::endl;
            std::cout<<"rotation: "<<std::endl<<der.first.rotation()<<std::endl<<std::endl;
        };
        
        numeric::Cluster3Geometry<K, TDirection3> clGeom;
        clGeom.init(sys);
        
        BOOST_CHECK(clGeom.parameters().size()==0);
        BOOST_CHECK(clGeom.derivatives().size()==0);
        
        clGeom.setBaseGeometry(&cluster);
        clGeom.recalculate();
        
        BOOST_CHECK(clGeom.parameters().size()==6);
        BOOST_CHECK(clGeom.derivatives().size()==6);
        
        cluster.addClusterGeometry(&clGeom);
        cluster.recalculate();
        
        //let's test the derivatives and see if we calculate them correct
        DerivativeTest::isCorrect(clGeom, 
                   std::bind(&numeric::Cluster3Geometry<K, TDirection3>::recalculate, &clGeom));
        
    }
    catch(boost::exception& x) {
        BOOST_FAIL(*boost::get_error_info<error_message>(x));
    }
    catch(std::exception& x) {
         BOOST_FAIL("Unknown exception");
    }
};


BOOST_AUTO_TEST_SUITE_END();
