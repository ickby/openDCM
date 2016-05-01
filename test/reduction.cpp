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

    //build up an example reduction tree
    dcm::symbolic::reduction::GeometryEdgeReductionTree<K, TDirection3, TDirection3> tree;
    
    //add a dependend geometry node 
    dcm::symbolic::reduction::Node node = tree.getTreeNode<PointLineGlider>();
    
    //try to add a edge
    tree.getSourceNode().connect(node, [](dcm::symbolic::reduction::TreeWalker* walker)->bool{
        
        auto cwalker = static_cast<dcm::symbolic::reduction::ConstraintWalker<K, TDirection3, TDirection3>*>(walker);
        
        
        return true;
    }
    );
}


BOOST_AUTO_TEST_SUITE_END();
