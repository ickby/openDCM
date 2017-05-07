/*
    openDCM, dimensional constraint manager
    Copyright (C) 2016  Stefan Troeger <stefantroeger@gmx.net>

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

#include <opendcm/core/constraint.hpp>
#include <opendcm/core/system.hpp>
#include <opendcm/core/reduction.hpp>

#include <boost/exception/get_error_info.hpp>

using namespace dcm;
      
template<typename Kernel>
struct TDirection3 : public geometry::Geometry<Kernel, numeric::Vector<Kernel, 3>> {

    using geometry::Geometry<Kernel, numeric::Vector<Kernel, 3>>::m_storage;
    //using geometry::Geometry<Kernel, geometry::storage::Vector<3>>::operator=;
   
    auto value() -> decltype(fusion::at_c<0>(m_storage)){
        return fusion::at_c<0>(m_storage);
    };
};

template<typename Kernel>
struct TScalar : public geometry::Geometry<Kernel, typename Kernel::Scalar> {

    using geometry::Geometry<Kernel, typename Kernel::Scalar>::m_storage;
   
    auto value() -> decltype(fusion::at_c<0>(m_storage)){
        return fusion::at_c<0>(m_storage);
    };
};

template<typename Kernel>
struct PointLineGlider : public numeric::DependendGeometry<Kernel, TDirection3, TScalar, double> {

    using Inherited = numeric::DependendGeometry<Kernel, TDirection3, TScalar, double>;

    CALCULATE() {
        Inherited::DependendGeometry::calculate();
        Inherited::value() = Inherited::input().value().norm()*scale();
    };
    
    auto scale() -> decltype(fusion::at_c<0>(Inherited::m_parameterStorage)){
        return fusion::at_c<0>(Inherited::m_parameterStorage);
    };
};

template<typename  Kernel>
struct FixedPoint : public numeric::DependendGeometry<Kernel, TDirection3, TDirection3> {
     
    using Inherited = numeric::DependendGeometry<Kernel, TDirection3, TDirection3>;

    CALCULATE() {
        Inherited::calculate();
        Inherited::output() = Inherited::input();
    };
};


namespace dcm { namespace numeric {
 
template<typename Kernel>
struct BinaryConstraint<Kernel, dcm::Distance, TDirection3<Kernel>, TDirection3<Kernel>> 
    : public BinaryConstraintBase<Kernel, dcm::Distance, TDirection3<Kernel>, TDirection3<Kernel>>  {

    typedef BinaryConstraintBase<Kernel, dcm::Distance, TDirection3<Kernel>, TDirection3<Kernel>>  Inherited;
    typedef typename Kernel::Scalar                 Scalar;
    typedef typename Inherited::Vector              Vector;
    typedef typename Inherited::Geometry1           Geometry1;
    typedef typename Inherited::Derivative1         Derivative1;
    typedef typename Inherited::Geometry2           Geometry2;
    typedef typename Inherited::Derivative2         Derivative2;
    
    BinaryConstraint() {
    };
    
    Scalar calculateError(const TDirection3<Kernel>& g1, const TDirection3<Kernel>& g2) {return 1;};
    Scalar calculateGradientFirst(const TDirection3<Kernel>& g1, const TDirection3<Kernel>& g2, const Derivative1& dg1) {return 1;};
    Scalar calculateGradientSecond(const TDirection3<Kernel>& g1, const TDirection3<Kernel>& g2, const Derivative2& dg2) {return 1;};
    Vector calculateGradientFirstComplete(const TDirection3<Kernel>& g1, const TDirection3<Kernel>& g2) {return Eigen::Vector2d(1,2);};
    Vector calculateGradientSecondComplete(const TDirection3<Kernel>& g1, const TDirection3<Kernel>& g2) {return Eigen::Vector2d(1,2);};
};

template<typename Kernel>
struct BinaryConstraint<Kernel, dcm::Angle, TDirection3<Kernel>, TDirection3<Kernel>> : public BinaryConstraintBase<Kernel, dcm::Angle, TDirection3<Kernel>, TDirection3<Kernel>>  {

    typedef BinaryConstraintBase<Kernel, dcm::Angle, TDirection3<Kernel>, TDirection3<Kernel>>  Inherited;
    typedef typename Kernel::Scalar                 Scalar;
    typedef typename Inherited::Vector              Vector;
    typedef typename Inherited::Geometry1           Geometry1;
    typedef typename Inherited::Derivative1         Derivative1;
    typedef typename Inherited::Geometry2           Geometry2;
    typedef typename Inherited::Derivative2         Derivative2;
    
    BinaryConstraint() {
    };
        
    Scalar calculateError(const TDirection3<Kernel>& g1, const TDirection3<Kernel>& g2) {return 1;};
    Scalar calculateGradientFirst(const TDirection3<Kernel>& g1, const TDirection3<Kernel>& g2, const Derivative1& dg1) {return 1;};
    Scalar calculateGradientSecond(const TDirection3<Kernel>& g1, const TDirection3<Kernel>& g2, const Derivative2& dg2) {return 1;};
    Vector calculateGradientFirstComplete(const TDirection3<Kernel>& g1, const TDirection3<Kernel>& g2) {return Eigen::Vector2d(1,2);};
    Vector calculateGradientSecondComplete(const TDirection3<Kernel>& g1, const TDirection3<Kernel>& g2) {return Eigen::Vector2d(1,2);};
};

template<typename Kernel>
struct UnaryConstraint<Kernel, dcm::Fix, TDirection3<Kernel>> 
    : public UnaryConstraintBase<Kernel, dcm::Fix, TDirection3<Kernel>>  {

    typedef UnaryConstraintBase<Kernel, dcm::Fix, TDirection3<Kernel>>  Inherited;
    typedef typename Kernel::Scalar                 Scalar;
    typedef typename Inherited::Vector              Vector;
    typedef typename Inherited::Geometry            Geometry;
    typedef typename Inherited::Derivative          Derivative;
    typedef typename Inherited::InputEquation       InputEquation;
    
    bool   applyToEquation(std::shared_ptr<InputEquation> eqn) {return false;};
    Scalar calculateError(const TDirection3<Kernel>& g1) {return 1;};
    Scalar calculateGradient(const TDirection3<Kernel>& g1, const Derivative& dg1) {return 1;};
    Vector calculateGradientComplete(const TDirection3<Kernel>& g1) {return Eigen::Vector2d(1,2);};
};

}}

struct ReductionModule {

    typedef boost::mpl::int_<1> ID;

    template<typename Final, typename Stacked>
    struct type : public Stacked {

        DCM_MODULE_ADD_GEOMETRIES(Stacked, (TDirection3)(TScalar));
    };
};

typedef dcm::Eigen3Kernel<double> K;
typedef dcm::System<K, ReductionModule> Sys;

template<typename Constraint, typename utilities::non_floating<typename Constraint::PimaryOptionType>::type Option>
using Equal = reduction::ConstraintEqualValue<K, Constraint, Option>;

BOOST_AUTO_TEST_SUITE(Reduction);

BOOST_AUTO_TEST_CASE(tree) {
        
    typedef dcm::reduction::GeometryEdgeReductionGraph<K, TDirection3, TDirection3> Tree;

    //build the data 
    TDirection3<K> g1, g2;
    auto sg1 = new dcm::symbolic::TypeGeometry<TDirection3<K>>();
    auto sg2 = new dcm::symbolic::TypeGeometry<TDirection3<K>>();
    sg1->setPrimitive(g1);
    sg2->setPrimitive(g2);
    
    auto c1 = new dcm::symbolic::TypeConstraint<dcm::Distance>();
    c1->setPrimitive(dcm::distance = 1.);
    
    auto c2 = new dcm::symbolic::TypeConstraint<dcm::Angle>();
    c2->setPrimitive(dcm::angle=2.);
    Tree::BinarySymbolicVector cvec;
    cvec.push_back(std::make_tuple(c1, sg1, sg2));
    cvec.push_back(std::make_tuple(c2, sg1, sg2));
    
    Tree::UnarySymbolicVector vvec;
    
    //build up an example reduction tree
    Tree tree;
    
    //add a dependend geometry node 
    std::shared_ptr<reduction::Node> node = tree.getTreeNode<PointLineGlider<K>>();
    
    //connect the node with a custom connection
    tree.sourceNode()->connect(node, 
        [](dcm::reduction::TreeWalker* walker)->symbolic::Constraint* {
            
            auto cwalker = static_cast<dcm::reduction::SourceTargetWalker<K, TDirection3, TDirection3>*>(walker);
            auto dist = cwalker->getConstraint<dcm::Distance>(dcm::Distance::index(),
                                                            dcm::Distance::Arity);
            auto val = dist->getPrimitive().distance();
            if(dist && std::abs((dist->getPrimitive().distance()-1))<1e-6) {
                    
                cwalker->acceptConstraint(dist);        
                return dist;
            }            
            return nullptr;
        }, 
        [](dcm::reduction::TreeWalker* walker, symbolic::Constraint*) {}
    );       
    
    //apply should execute both nodes and the connection
    auto walker = static_cast<dcm::reduction::SourceTargetWalker<K, TDirection3, TDirection3>*>(tree.apply(sg1, sg2, cvec, vvec));
    
    BOOST_CHECK(!walker->getInputEquation());
    BOOST_CHECK(walker->getInitialNode() == tree.sourceNode());
    BOOST_CHECK(walker->getFinalNode()  == tree.getTreeNode<PointLineGlider<K>>());
    BOOST_CHECK_EQUAL(walker->binaryConstraintPool().size(), 1);
    
    //lets test conditional constraints 
    delete walker;
    c1->getPrimitive().distance() = 2;
    auto fixed = tree.getTreeNode<FixedPoint<K>>();

    tree.sourceNode()->connectConditional<reduction::ConstraintEqualValue<K, dcm::Angle, 2>>(fixed, 
                                                                [](reduction::ConstraintWalker<K>* w, const dcm::Angle& angle){}); 
    
    walker = static_cast<dcm::reduction::SourceTargetWalker<K, TDirection3, TDirection3>*>(tree.apply(sg1, sg2, cvec, vvec));

    BOOST_CHECK(!walker->getInputEquation());
    BOOST_CHECK(walker->getInitialNode() == tree.sourceNode());
    BOOST_CHECK(walker->getFinalNode()  == tree.getTreeNode<FixedPoint<K>>());
    BOOST_CHECK_EQUAL(walker->binaryConstraintPool().size(), 1);
}

BOOST_AUTO_TEST_CASE(convertion) {
    
    try{
    //build up a small test graph
    typedef typename Sys::Graph Graph;
    std::shared_ptr<Graph> g = std::shared_ptr<Graph>(new Graph);
    auto v1 = g->addVertex();
    auto v2 = g->addVertex();
    auto e1  = g->addEdge(fusion::at_c<0>(v1), fusion::at_c<0>(v2));
    auto e2  = g->addEdge(fusion::at_c<0>(v1), fusion::at_c<0>(v2));
    
    //build the data 
    TDirection3<K> g1, g2;
    g1.value() << 1,2,3;
    g2.value() << 4,5,6;
    auto sg1 = new dcm::symbolic::TypeGeometry<TDirection3<K>>();
    auto sg2 = new dcm::symbolic::TypeGeometry<TDirection3<K>>();
    sg1->setPrimitive(g1);
    sg2->setPrimitive(g2);
    g->setProperty<symbolic::GeometryProperty>(fusion::at_c<0>(v1), sg1);
    g->setProperty<symbolic::GeometryProperty>(fusion::at_c<0>(v2), sg2);

    auto c1 = new dcm::symbolic::TypeConstraint<dcm::Distance>();
    c1->setPrimitive(dcm::distance);
    auto c2 = new dcm::symbolic::TypeConstraint<dcm::Angle>();
    c2->setPrimitive(dcm::angle);
    g->setProperty<symbolic::ConstraintProperty>(fusion::at_c<1>(e1), c1);
    g->setProperty<symbolic::ConstraintProperty>(fusion::at_c<1>(e2), c2);

    symbolic::NumericConverter<K, typename Sys::GeometryList, typename Sys::BinaryConstraintList,
                               typename Sys::UnaryConstraintList, typename Sys::Graph> reducer;
                               
    reducer.setupEquationHandler(g, fusion::at_c<0>(e1));
    
    //check the result
    numeric::EquationHandler<K>* builder = g->getProperty<numeric::EquationHandlerProperty<K>>(fusion::at_c<0>(e1));
    BOOST_CHECK(builder != nullptr);
    
    auto eg1 = std::static_pointer_cast<numeric::Equation<K, TDirection3<K>>>(builder->createGeometry(fusion::at_c<0>(v1)));
    auto eg2 = std::static_pointer_cast<numeric::Equation<K, TDirection3<K>>>(builder->createGeometry(fusion::at_c<0>(v2)));
    BOOST_CHECK( eg1->output().value().isApprox(Eigen::Vector3d(1,2,3), 1e-9) );
    BOOST_CHECK( eg2->output().value().isApprox(Eigen::Vector3d(4,5,6), 1e-9) );
    
    auto ecv = builder->createBinaryEquations(eg1, eg2);
    //BOOST_CHECK_EQUAL(ecv.size(),2);
    auto ec = std::static_pointer_cast<numeric::BinaryEquation<K, TDirection3<K>, TDirection3<K>, double>>(ecv.front());
    BOOST_REQUIRE( ec );
    BOOST_CHECK( ec->firstInput().value().isApprox(Eigen::Vector3d(1,2,3), 1e-9) );
    BOOST_CHECK( ec->secondInput().value().isApprox(Eigen::Vector3d(4,5,6), 1e-9) );
    
    }
    catch(const boost::exception& e) {
        std::cout<<"unexpected error " << *boost::get_error_info<boost::errinfo_errno>(e)
                    << " raised: " << boost::get_error_info<dcm::error_message>(e)->c_str()<<std::endl;
    };
}

BOOST_AUTO_TEST_CASE(unary) {
 
    typedef dcm::reduction::GeometryEdgeReductionGraph<K, TDirection3, TDirection3> Tree;

    //build the data 
    TDirection3<K> g1, g2;
    auto sg1 = new dcm::symbolic::TypeGeometry<TDirection3<K>>();
    auto sg2 = new dcm::symbolic::TypeGeometry<TDirection3<K>>();
    sg1->setPrimitive(g1);
    sg2->setPrimitive(g2);
    
    auto c1 = new dcm::symbolic::TypeConstraint<dcm::Distance>();
    c1->setPrimitive(dcm::distance);
    c1->getPrimitive().distance() = 0;
    
    auto c2 = new dcm::symbolic::TypeConstraint<dcm::Fix>();
    c2->setPrimitive(dcm::fix=dcm::Fixables::pointY);
    c2->getPrimitive().fixed() = dcm::Fixables::pointY;
    Tree::BinarySymbolicVector cvec;
    cvec.push_back(std::make_tuple(c1, sg1, sg2));
   
    Tree::UnarySymbolicVector vvec;
    vvec.push_back(std::make_tuple(c2, sg2));
    
    //build up an example reduction tree
    Tree tree;
    
    //add a dependend geometry node 
    std::shared_ptr<reduction::Node> node = tree.getTreeNode<PointLineGlider<K>>();
    
    //connect the node with a custom connection
    tree.sourceNode()->connect(node, 
        [](dcm::reduction::TreeWalker* walker)-> symbolic::Constraint* {
            
            auto cwalker = static_cast<dcm::reduction::SourceTargetWalker<K, TDirection3, TDirection3>*>(walker);
            auto dist = cwalker->getConstraint<dcm::Fix>(dcm::Fix::index(),
                                                        dcm::Fix::Arity);
            if(dist && dist->getPrimitive().fixed()==dcm::Fixables::pointY) {
                    
                cwalker->acceptConstraint(dist);
                return dist;
            }            
            return nullptr;
        },
        [](dcm::reduction::TreeWalker* walker, symbolic::Constraint*) {}
    );       
    
    //apply should execute both nodes and the connection
    auto walker = static_cast<dcm::reduction::SourceTargetWalker<K, TDirection3, TDirection3>*>(tree.apply(sg1, sg2, cvec, vvec));
    
    BOOST_CHECK(!walker->getInputEquation());
    BOOST_CHECK(walker->getInitialNode() == tree.sourceNode());
    BOOST_CHECK(walker->getFinalNode()  == tree.getTreeNode<PointLineGlider<K>>());
    BOOST_CHECK_EQUAL(walker->binaryConstraintPool().size(), 1);
    BOOST_CHECK_EQUAL(walker->unaryConstraintPool().size(), 0);    
}

BOOST_AUTO_TEST_SUITE_END();
