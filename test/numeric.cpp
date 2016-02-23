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

#include "opendcm/core/equations.hpp"
#include "opendcm/core/geometry.hpp"
#include "opendcm/core/kernel.hpp"

#include <boost/fusion/include/at.hpp>

using namespace dcm;

typedef dcm::Eigen3Kernel<double> K;

//a default equation wich uses own parameters
struct TestEquation : public numeric::Equation<K, Eigen::Vector3d> {
    
    typedef numeric::Equation<K, Eigen::Vector3d> Base;
    
    TestEquation() : Base(), m_map(nullptr) {};
    TestEquation(const Eigen::Vector3d& val) : Base(val), m_map(nullptr) {};
    
    virtual void init(numeric::LinearSystem< K >& sys) {
        
        Base::m_parameters = sys.mapParameter(m_map);
        
        for(auto& param : Base::m_parameters) 
            Base::m_derivatives.push_back(std::make_pair(Eigen::Vector3d(), param));
    };
    
    CALCULATE() {
    
        //nothing to calculate, only updating
        output() = m_map;
        
        //derivatives are simple. should actually already be defined in init, but this is test only
        for(int i=0; i<3; ++i)
            Base::m_derivatives[i].first(i) = 1;
        
        //this equation is not used as error function, hence no derivtives and result update
    };
    
private:
    Eigen::Map<Eigen::Vector3d> m_map;
};


template<typename Kernel, bool MappedType = true>
struct TDirection3 : public geometry::Geometry<Kernel, MappedType,
            geometry::storage::Vector<3>> {

    using geometry::Geometry<Kernel, MappedType, geometry::storage::Vector<3>>::m_storage;
    using geometry::Geometry<Kernel, MappedType, geometry::storage::Vector<3>>::operator=;
   
    auto value() -> decltype(fusion::at_c<0>(m_storage)){
        return fusion::at_c<0>(m_storage);
    };
};
 
template<typename Kernel, bool MappedType = true>
struct TMatrix3 : public geometry::Geometry<Kernel, MappedType,
            geometry::storage::Matrix<3,3>> {

    using geometry::Geometry<Kernel, MappedType, geometry::storage::Matrix<3,3>>::m_storage;
    
    auto value() -> decltype(fusion::at_c<0>(m_storage)){
        return fusion::at_c<0>(m_storage);
    };
};

template<typename Kernel, bool MappedType = true>
struct TCylinder3 : public geometry::Geometry<Kernel, MappedType,
            geometry::storage::Vector<3>, geometry::storage::Parameter, geometry::storage::Vector<3>> {

    typedef geometry::Geometry<Kernel, MappedType, geometry::storage::Vector<3>, 
                geometry::storage::Parameter, geometry::storage::Vector<3>> Inherited;
    using Inherited::m_storage;
    
    auto point() -> decltype(fusion::at_c<0>(m_storage)){
        return fusion::at_c<0>(m_storage);
    };
    
    auto direction() -> decltype(fusion::at_c<2>(m_storage)){
        return fusion::at_c<2>(m_storage);
    };
    
    typename Kernel::Scalar& radius() {
        return Inherited::rmPtr(fusion::at_c<1>(m_storage));
    };
};

using namespace dcm;



BOOST_AUTO_TEST_SUITE(Numeric_test_suit);

BOOST_AUTO_TEST_CASE(equations) {
    
    //see if we can construct it
    numeric::Equation<K, double> eqn1(2), eqn2(2.5);
    numeric::Equation<K, Eigen::Vector3d> eqn3, eqn4(Eigen::Vector3d(1,2,3));
    
    //accessing?
    double d1 = eqn2.output();
    double d2 = eqn2;
    Eigen::Vector3d v1 = eqn4.output();
    Eigen::Vector3d v2 = eqn3;
    
    //using special functions of the return type?
    double d3 = eqn4.norm();
    
    BOOST_CHECK_CLOSE(d1, 2.5, 1e-9);
    BOOST_CHECK_CLOSE(d2, 2.5, 1e-9);
    BOOST_CHECK_CLOSE(v1(0), 1, 1e-9);
    BOOST_CHECK_CLOSE(v1(1), 2, 1e-9);
    BOOST_CHECK_CLOSE(v1(2), 3, 1e-9);   
    
    
    //check unary and binary equations
    numeric::LinearSystem<K> sys(20,20); 
    
    auto d = std::make_shared<numeric::Equation<K, double>>(eqn1);
    auto v = std::shared_ptr<TestEquation>(new TestEquation());
    
    auto d_v = numeric::makeUnaryEquation<K, double, Eigen::Vector3d>  (
        [](const double& d, Eigen::Vector3d& v) {
            v << d, 2*d, 3*d;    
        },        
        [](const double& d, const double& dd, Eigen::Vector3d& v) {
            v << 1*dd,2*dd,3*dd;
        }
    );
    
    auto vv_d = numeric::makeBinaryEquation<K, Eigen::Vector3d, Eigen::Vector3d, double>  (
        [](const Eigen::Vector3d& d1, const Eigen::Vector3d& d2, double& v) {
            v = d1.transpose() * d2;    
        },        
        [](const Eigen::Vector3d& d1, const Eigen::Vector3d& d2, const Eigen::Vector3d& dd1, double& v) {
            v = dd1.transpose() * d2; 
        },
        [](const Eigen::Vector3d& d1, const Eigen::Vector3d& d2, const Eigen::Vector3d& dd2, double& v) {
            v = d1.transpose() * dd2; 
        }
    );
    
    d_v->prepend(d);
    d_v->execute();
    BOOST_CHECK_SMALL((*d_v-Eigen::Vector3d(2,4,6)).norm(), 1e-9);
    BOOST_CHECK(d_v->parameters().size() == 0);
    BOOST_CHECK(d_v->derivatives().size() == 0);
    
    d_v->init(sys);    
    BOOST_CHECK(d_v->parameters().size() == 0);
    BOOST_CHECK(d_v->derivatives().size() == 0);
    
    vv_d->prepend(v, d_v);
    BOOST_CHECK(vv_d->parameters().size() == 0);
    BOOST_CHECK(vv_d->derivatives().size() == 0);
    
    vv_d->init(sys);
    BOOST_CHECK_EQUAL(vv_d->parameters().size(), 0);
    BOOST_CHECK_EQUAL(vv_d->derivatives().size(),3);
}


BOOST_AUTO_TEST_CASE(geometry) {
  
    TDirection3<K>        basic; //cannot use it as it is a nullptr initialized map
    TDirection3<K, false> basic_vec;
    
    basic_vec.value()[0] = 5;    
    BOOST_CHECK(basic_vec.value()[0] == 5);
    
    basic_vec.value() = Eigen::Vector3d(1,2,3);    
    BOOST_CHECK(basic_vec.value() == Eigen::Vector3d(1,2,3));
    
    numeric::LinearSystem<K> sys(20,20);    
    numeric::Geometry<K, TDirection3> dirGeom;
    
    BOOST_CHECK(dirGeom.newParameterCount() == 3);
    
    dirGeom.init(sys);
    BOOST_REQUIRE(dirGeom.parameters().size() == 3);
    BOOST_REQUIRE(dirGeom.derivatives().size() == 3);
    
    typedef numeric::Geometry<K, TDirection3>::DerivativePack    TDir3Derivative;
    typedef numeric::Geometry<K, TDirection3>::ParameterIterator TDir3ParIt;

    int c = 0;
    TDir3ParIt it = dirGeom.parameters().begin();
    for(TDir3Derivative& der : dirGeom.derivatives()) {
        BOOST_CHECK(der.first.value()(c)==1);
        BOOST_CHECK(der.first.value().sum()==1);
        BOOST_CHECK(der.second == *(it++));
        BOOST_CHECK(der.second.Value != nullptr);
        BOOST_CHECK(der.second.Value != 0);
        ++c;
    };
    
    numeric::Geometry<K, TCylinder3> cylGeom;
    
    cylGeom.init(sys);
    BOOST_REQUIRE(cylGeom.parameters().size() == 7);
    BOOST_REQUIRE(cylGeom.derivatives().size() == 7);
   
    typedef numeric::Geometry<K, TCylinder3>::DerivativePack        TCyl3Derivative;
    typedef numeric::Geometry<K, TCylinder3>::ParameterIterator  TCyl3ParIt;
    
    c=0;
    Eigen::VectorXd res(7);
    TCyl3ParIt cylIt = cylGeom.parameters().begin();
    for(TCyl3Derivative& der : cylGeom.derivatives()) {
        res.setZero();
        res(c) = 1;
        
        BOOST_CHECK(der.first.point()(0)==res(0));
        BOOST_CHECK(der.first.point()(1)==res(1));
        BOOST_CHECK(der.first.point()(2)==res(2));
        BOOST_CHECK(der.first.radius()==res(3));
        BOOST_CHECK(der.first.direction()(0)==res(4));
        BOOST_CHECK(der.first.direction()(1)==res(5));
        BOOST_CHECK(der.first.direction()(2)==res(6));
        BOOST_CHECK(der.second == *(cylIt++));
        BOOST_CHECK(der.second.Value != nullptr);
        BOOST_CHECK(der.second.Value != 0);
        
        c++;
    };
    
    //let's see if the mapping works
    sys.parameter().head<3>() = Eigen::Vector3d(7.1,8.2,9.3);
    dirGeom.execute();
    BOOST_CHECK(dirGeom.value().isApprox(Eigen::Vector3d(7.1,8.2,9.3)));
    
    sys.parameter()(6) = 3.3;
    cylGeom.execute();
    BOOST_CHECK(cylGeom.radius() == (3.3));
    
    sys.parameter().segment<3>(7) = Eigen::Vector3d(5.5,6.6,7.7);
    cylGeom.execute();
    BOOST_CHECK(cylGeom.direction().isApprox(Eigen::Vector3d(5.5,6.6,7.7)));
    
    //see if anything was overriden
    BOOST_CHECK(dirGeom.value().isApprox(Eigen::Vector3d(7.1,8.2,9.3)));
    BOOST_CHECK(cylGeom.radius() == (3.3));
    BOOST_CHECK(cylGeom.direction().isApprox(Eigen::Vector3d(5.5,6.6,7.7)));
    BOOST_CHECK(cylGeom.point().isApprox(Eigen::Vector3d(0,0,0)));
    
    
    //and check matrix types
    numeric::Geometry<K, TMatrix3> matGeom;
    matGeom.init(sys);
    
    matGeom.value() << 1,2,3,4,5,6,7,8,9;
    BOOST_CHECK(matGeom.value().col(0).isApprox(Eigen::Vector3d(1,4,7)));
    BOOST_CHECK(matGeom.value().col(1).isApprox(Eigen::Vector3d(2,5,8)));
    BOOST_CHECK(matGeom.value().col(2).isApprox(Eigen::Vector3d(3,6,9)));
    
    //check assign and copy-constructability
    TDirection3<K, false> basic_vec_2(basic_vec);
    BOOST_CHECK(basic_vec_2.value().isApprox(basic_vec.value()));
    basic_vec.value() << 5,4,3;
    basic_vec_2 = basic_vec;
    BOOST_CHECK(basic_vec_2.value().isApprox(basic_vec.value()));
};
    

BOOST_AUTO_TEST_CASE(parameter_geometry) {

    numeric::LinearSystem<K> sys(10,10); 
    Eigen::VectorXd init = sys.parameter();
    
    numeric::ParameterGeometry<K, TCylinder3, dcm::geometry::storage::Parameter> cylGeom;
    
    BOOST_CHECK(cylGeom.newParameterCount() == 1);
    
    cylGeom.init(sys);       
    
    //check default constructed derivatives
    BOOST_CHECK(cylGeom.parameters().size()==1);
    BOOST_CHECK(cylGeom.derivatives().size()==1);
    
    //see if the mapping to the internal storage worked
    cylGeom.point() = Eigen::Vector3d(1,2,3);    
    BOOST_CHECK(cylGeom.point().isApprox(Eigen::Vector3d(1,2,3)));
    BOOST_CHECK(sys.parameter().isApprox(init));
       
    numeric::ParameterGeometry<K, TCylinder3, dcm::geometry::storage::Parameter, 
                    dcm::geometry::storage::Vector<3>> cyl2Geom;
                   
    BOOST_CHECK(cyl2Geom.newParameterCount() == 4);
                    
    cyl2Geom.init(sys);
    BOOST_CHECK(cyl2Geom.newParameterCount() == 4);
    BOOST_CHECK(cyl2Geom.parameters().size()==4);
    BOOST_CHECK(cyl2Geom.derivatives().size()==4);
    
    //check the polymorphism
    numeric::Equation<K, TCylinder3<K, false>>* Geom = &cylGeom;
    
    BOOST_CHECK(Geom->point().isApprox(Eigen::Vector3d(1,2,3)));
    BOOST_CHECK(sys.parameter().isApprox(init));
    
    cylGeom.direction() = Eigen::Vector3d(4,5,6);
    BOOST_CHECK(Geom->direction().isApprox(Eigen::Vector3d(4,5,6)));
    BOOST_CHECK(sys.parameter().isApprox(init));
    
    Geom->radius() = 7;
    BOOST_CHECK(cylGeom.radius() == 7);
    BOOST_CHECK(sys.parameter().isApprox(init));
    
    //check if the parameter mapping worked
    BOOST_CHECK(cylGeom.parameters()[0] == cylGeom.derivatives()[0].second);
    *cylGeom.parameters()[0].Value = 5;
    BOOST_CHECK(cylGeom.derivatives()[0].second.Value != nullptr);
    BOOST_CHECK(*cylGeom.derivatives()[0].second.Value == 5);
    
    //check assign and copy-constructability
    numeric::ParameterGeometry<K, TCylinder3, 
        dcm::geometry::storage::Parameter>  cylGeom2(cylGeom);
    BOOST_CHECK(cylGeom2.point().isApprox(cylGeom.point()));
    BOOST_CHECK(cylGeom2.direction().isApprox(cylGeom.direction()));
    BOOST_CHECK(cylGeom2.radius() == cylGeom.radius());
    cylGeom.point() << 0,0,3;
    cylGeom.direction() << -3,6,1;
    cylGeom.radius() = -7;
    cylGeom2 = cylGeom;
    BOOST_CHECK(cylGeom2.point().isApprox(cylGeom.point()));
    BOOST_CHECK(cylGeom2.direction().isApprox(cylGeom.direction()));
    BOOST_CHECK(cylGeom2.radius() = cylGeom.radius());
};


BOOST_AUTO_TEST_SUITE_END();
