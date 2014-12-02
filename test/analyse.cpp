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

#include <opendcm/core/geometry.hpp>
#include <opendcm/core/analyse.hpp>
#include <opendcm/core/system.hpp>
#include <boost/multi_array.hpp>

using namespace dcm;

typedef Eigen3Kernel<double> K;

template<typename Kernel = K, bool Map = false>
struct TDirection3 : public geometry::Geometry<Kernel, Map,
        geometry::storage::Vector<3>> {

    typedef typename Kernel::Scalar Scalar;
    using geometry::Geometry<Kernel, Map, geometry::storage::Vector<3>>::m_storage;

    auto value() -> decltype(fusion::at_c<0>(m_storage)) {
        return fusion::at_c<0>(m_storage);
    };

    TDirection3<Kernel, Map>& transform(const Eigen::Transform<Scalar, 3, Eigen::AffineCompact>& t) {
        value() = t.rotation()*value();
        return *this;
    };

    TDirection3<Kernel, Map>  transformed(const Eigen::Transform<Scalar, 3, Eigen::AffineCompact>& t) {
        TDirection3<Kernel, Map> copy(*this);
        copy.transform(t);
        return copy;
    };
};

struct TestModule1 {
    typedef boost::mpl::int_<1> ID;

    template<typename Final, typename Stacked>
    struct type : public Stacked {
        DCM_MODULE_ADD_GEOMETRIES(Stacked, (TDirection3) )
        DCM_MODULE_ADD_VERTEX_PROPERTIES(Stacked, (symbolic::GeometryProperty))
    };
};

typedef System<TestModule1> TestSystem;

struct value {
    typedef int type;
};

struct symbol {
    typedef int type;
};

typedef graph::ClusterGraph<mpl::vector1<value>, mpl::vector1<symbol>,
        mpl::vector0<>, mpl::vector0<> > Graph;

        
struct Node1 : public symbolic::GeometryNode<K, TDirection3, TDirection3> {
    
    virtual void execute(symbolic::GeometryTreeWalker< K, TDirection3, TDirection3 >* walker) {
        std::cout<<"excute"<<std::endl;
    };
};        

BOOST_AUTO_TEST_SUITE(analyse);

BOOST_AUTO_TEST_CASE(analyse_basic) {

    boost::multi_array<symbolic::EdgeReductionTree<TestSystem>*,2> reduction;
    reduction.resize(boost::extents[1][1]);
    reduction[0][0] = new symbolic::GeometryEdgeReductionTree<TestSystem, TDirection3, TDirection3>();
};

//     //build up the graph
//     boost::shared_ptr< Graph > g = boost::shared_ptr< Graph >(new Graph);
//     
//     dcm::GlobalVertex v1 = fusion::at_c<1>(g->addVertex());
//     dcm::GlobalVertex v2 = fusion::at_c<1>(g->addVertex());
//     
//     dcm::GlobalEdge e1 = fusion::at_c<1>(g->addEdge(v1,v2));
//     dcm::GlobalEdge e2 = fusion::at_c<1>(g->addEdge(v1,v2));
//     dcm::GlobalEdge e3 = fusion::at_c<1>(g->addEdge(v1,v2));
//     dcm::GlobalEdge e4 = fusion::at_c<1>(g->addEdge(v1,v2));
//     
//     g->setProperty<symbol>(e1, 1);
//     g->setProperty<symbol>(e2, 2);
//     g->setProperty<symbol>(e3, 1);
//     g->setProperty<symbol>(e4, 2);
//     
//     //build up the analyser tree

BOOST_AUTO_TEST_SUITE_END();
