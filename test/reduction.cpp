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

#include <opendcm/core/reduction.hpp>

using namespace dcm;

typedef dcm::Eigen3Kernel<double> K;
typedef dcm::graph::ClusterGraph<mpl::vector0<>, mpl::vector1<dcm::symbolic::ConstraintProperty>,
        mpl::vector1<dcm::symbolic::GeometryProperty>, mpl::vector0<> > Graph;
        
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


struct PointLineGlider : public numeric::DependendGeometry<K, TDirection3, TScalar, double> {
     
    CALCULATE() {
        DependendGeometry::calculate();
        value() = input().value().norm()*scale();
    };
    
    auto scale() -> decltype(fusion::at_c<0>(m_parameterStorage)){
        return fusion::at_c<0>(m_parameterStorage);
    };
};

BOOST_AUTO_TEST_SUITE(Reduction);

BOOST_AUTO_TEST_CASE(tree) {
    
    //build up a small test graph
    std::shared_ptr<Graph> g = std::shared_ptr<Graph>(new Graph);
    auto v1 = g->addVertex();
    auto v2 = g->addVertex();
    auto e  = g->addEdge(fusion::at_c<0>(v1), fusion::at_c<0>(v2));
    
    //build the data 
    TDirection3<K> g1, g2;
    auto sg1 = new dcm::symbolic::TypeGeometry<K, TDirection3>();
    auto sg2 = new dcm::symbolic::TypeGeometry<K, TDirection3>();
    sg1->setPrimitiveGeometry(g1);
    sg2->setPrimitiveGeometry(g2);
    auto c1 = new dcm::symbolic::TypeConstraint<dcm::Distance>();
    c1->setPrimitiveConstraint(dcm::distance);
    c1->setConstraintID(0);
    auto c2 = new dcm::symbolic::TypeConstraint<dcm::Angle>();
    c2->setPrimitiveConstraint(dcm::angle);
    c2->setConstraintID(1);
    std::vector<symbolic::Constraint*> cvec;
    cvec.push_back(c1);
    cvec.push_back(c2);
    
    //build up an example reduction tree
    dcm::symbolic::reduction::GeometryEdgeReductionTree<K, TDirection3, TDirection3> tree;
    
    //add a dependend geometry node 
    dcm::symbolic::reduction::Node* node = tree.getTreeNode<PointLineGlider>();
    
    //connect the node with a custom edge
    tree.getSourceNode()->connect(node, [](dcm::symbolic::reduction::TreeWalker* walker)->bool{
        
        auto cwalker = static_cast<dcm::symbolic::reduction::ConstraintWalker<K, TDirection3, TDirection3>*>(walker);
        auto dist = cwalker->getConstraint<dcm::Distance>(1);
        if(dist && dist->getPrimitiveConstraint().distance()==0) {
                
            cwalker->acceptConstraint(dist);        
            return true;
        }            
        return false;
    }
    );       
    
    //apply should execute both nodes and the edge
    auto walker = static_cast<dcm::symbolic::reduction::ConstraintWalker<K, TDirection3, TDirection3>*>(tree.apply(sg1, sg2, cvec));
    
    BOOST_CHECK(!walker->getCummulativeInputEquation());
    BOOST_CHECK(walker->getGeometry());
    BOOST_CHECK(walker->size() == 1);
    
}


BOOST_AUTO_TEST_SUITE_END();
