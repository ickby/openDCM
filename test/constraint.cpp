/*
    openDCM, dimensional constraint manager
    Copyright (C) 2014  Stefan Troeger <stefantroeger@gmx.net>

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
#include <boost/concept_check.hpp>

#include "opendcm/core/constraint.hpp"
#include <Eigen/Core>

typedef dcm::Eigen3Kernel<double> K;

//two vectors perpendicular, maybe the easiest constraints of them all
struct test_constraint1 : public dcm::constraint::Constraint<test_constraint1, 2, int> {
    using Constraint::operator=;
    test_constraint1(const int& i) : Constraint(i) {};
    
    int& radius() {return fusion::at_c<0>(m_storage);};
};

struct test_constraint2 : public dcm::constraint::Constraint<test_constraint2, 2, double, char> {
    using Constraint::operator=;
    test_constraint2(const double& d, const char& c) : Constraint(d, c) {};
    
    double& radius()    {return fusion::at_c<0>(m_storage);};
    char&   direction() {return fusion::at_c<1>(m_storage);};
};

template<typename Kernel>
struct TPoint3 : public dcm::geometry::Geometry<Kernel, dcm::numeric::Vector<Kernel, 3>> {

    using dcm::geometry::Geometry<Kernel, dcm::numeric::Vector<Kernel, 3>>::m_storage;
   
    auto value() -> decltype(fusion::at_c<0>(m_storage)){
        return fusion::at_c<0>(m_storage);
    };
};

namespace dcm {
namespace numeric {

template<typename Kernel>
struct BinaryConstraint<Kernel, dcm::Distance, TPoint3<Kernel>, TPoint3<Kernel>> : public BinaryConstraintBase<Kernel, dcm::Distance, TPoint3<Kernel>, TPoint3<Kernel>> {
  
    typedef BinaryConstraintBase<Kernel, dcm::Distance, TPoint3<Kernel>, TPoint3<Kernel>>  Inherited;
    typedef typename Kernel::Scalar                 Scalar;
    typedef typename Inherited::Vector              Vector;
    typedef typename Inherited::Geometry1           Geometry1;
    typedef typename Inherited::Derivative1         Derivative1;
    typedef typename Inherited::Geometry2           Geometry2;
    typedef typename Inherited::Derivative2         Derivative2;
    
    BinaryConstraint() {
    };
    
    Scalar calculateError(Geometry1& g1, Geometry2& g2) {        
        return (g1.value()-g2.value()).norm() - Inherited::distance();
    };

    Scalar calculateGradientFirst(Geometry1& g1, Geometry2& g2, Derivative1& dg1) {        
        return (g1.value()-g2.value()).dot(dg1.value()) / (g1.value()-g2.value()).norm();
    };

    Scalar calculateGradientSecond(Geometry1& g1, Geometry2& g2, Derivative2& dg2) {        
        return (g1.value()-g2.value()).dot(-dg2.value()) / (g1.value()-g2.value()).norm();
    };

    Vector calculateGradientFirstComplete(Geometry1& g1, Geometry2& g2) {
        return (g1.value()-g2.value()) / (g1.value()-g2.value()).norm();
    };

    Vector calculateGradientSecondComplete(Geometry1& g1, Geometry2& g2) {
        return (g1.value()-g2.value()) / (g1.value()-g2.value()).norm();
    };
};

} //numeric
} //dcm


template<typename T>
void pretty(T t) {
    std::cout << __PRETTY_FUNCTION__ << std::endl;
};

BOOST_AUTO_TEST_SUITE(Constraint_test_suit);

BOOST_AUTO_TEST_CASE(primitive) {
    
    test_constraint1 test1(3);
    test_constraint2 test2(0.1, 'b');
    test1 = 2;
    test2(0.2,'a');
    
    BOOST_CHECK(test1.radius() == 2);
    BOOST_CHECK(test2.radius() == 0.2);
    BOOST_CHECK(test2.direction() == 'a');
    
    test1.setDefault();
    test2.setDefault();
    
    BOOST_CHECK(test1.radius() == 3);
    BOOST_CHECK(test2.radius() == 0.1);
    BOOST_CHECK(test2.direction() == 'b');

    test1(test_constraint1(2));
    test2(test_constraint2(0.2,'a'));
     
    BOOST_CHECK(test1.radius() == 2);
    BOOST_CHECK(test2.radius() == 0.2);
    BOOST_CHECK(test2.direction() == 'a');
    
    test1 = test_constraint1(3);
    test2 = test_constraint2(0.1, 'b');
    
    BOOST_CHECK(test1.radius() == 3);
    BOOST_CHECK(test2.radius() == 0.1);
    BOOST_CHECK(test2.direction() == 'b');    
}

BOOST_AUTO_TEST_CASE(numeric) {
    
   dcm::numeric::LinearSystem<K> sys(20,20);  
   std::shared_ptr<dcm::numeric::Geometry<K, TPoint3>> p1(new dcm::numeric::Geometry<K, TPoint3>());
   std::shared_ptr<dcm::numeric::Geometry<K, TPoint3>> p2(new dcm::numeric::Geometry<K, TPoint3>);

   p1->init(sys);
   p2->init(sys);   
   p1->value() = Eigen::Vector3d(1,0,0);
   p2->value() = Eigen::Vector3d(0,0,0);
   
   typedef dcm::numeric::BinaryConstraintSimplifiedEquation<K, dcm::Distance, TPoint3<K>, TPoint3<K>>        ggc;
   typedef dcm::numeric::BinaryConstraintComplexEquation<K, dcm::Distance, TPoint3<K>, TPoint3<K>>           ccc;
   typedef dcm::numeric::BinaryConstraintSimplifiedComplexEquation<K, dcm::Distance, TPoint3<K>, TPoint3<K>> gcc;
   typedef dcm::numeric::BinaryConstraintComplexSimplifiedEquation<K, dcm::Distance, TPoint3<K>, TPoint3<K>> cgc;
   
   std::shared_ptr<ggc> gg_constraint(new ggc());
   std::shared_ptr<ccc> cc_constraint(new ccc());
   std::shared_ptr<gcc> gc_constraint(new gcc());
   std::shared_ptr<cgc> cg_constraint(new cgc());

   gg_constraint->setInputEquations(p1, p2);
   gg_constraint->init(sys);
   cc_constraint->setInputEquations(p1, p2);
   cc_constraint->init(sys);
   cg_constraint->setInputEquations(p1, p2);
   cg_constraint->init(sys);
   gc_constraint->setInputEquations(p1, p2);
   gc_constraint->init(sys);
   
   BOOST_CHECK_NO_THROW(gg_constraint->calculate());
   BOOST_CHECK_NO_THROW(cg_constraint->operator()());
   BOOST_CHECK_NO_THROW(cc_constraint->execute());
   BOOST_CHECK_NO_THROW(gc_constraint->calculate());
   
   BOOST_CHECK(gg_constraint->getResidual() == 1);

}

BOOST_AUTO_TEST_SUITE_END();
