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

#ifndef DCM_ANALYSE_H
#define DCM_ANALYSE_H

#include <Eigen/Core>

#include "geometry.hpp"
#include "clustergraph.hpp"

namespace dcm {
namespace symbolic {
    
struct ReductionResult {
  
    virtual void getGeometry() {};
};
   
struct ResultProperty {
    typedef ReductionResult* type;
    struct default_value {
        ReductionResult* operator()() {
            return nullptr;
        };
    };
};

template<typename Kernel, template<class, bool> class G1, template<class, bool> class G2>
struct Node;

template<typename Kernel, template<class, bool> class G1, template<class, bool> class G2>
struct GeometryTreeWalker;
    
template<typename Kernel, template<class, bool> class G1, template<class, bool> class G2>
struct GeometryEdge {

    Node<Kernel, G1, G2>* start;
    Node<Kernel, G1, G2>* end;
    
    virtual bool apply(symbolic::GeometryTreeWalker<Kernel, G1, G2>* walker) const = 0;
};
   
template<typename Kernel, template<class, bool> class G1, template<class, bool> class G2>
struct GeometryNode  {
    
    virtual void execute(symbolic::GeometryTreeWalker<Kernel, G1, G2>* walker) const = 0;
    
protected:    
    bool applyWalker(GeometryTreeWalker<Kernel, G1, G2>* walker) {
        for(const GeometryEdge<Kernel, G1, G2>& e : edges) {
            if(e.apply(walker))
                break;
        };
        return true;
    };
    
    std::vector<GeometryEdge<Kernel, G1, G2>> edges;
};

template<typename Kernel, template<class, bool> class G1, template<class, bool> class G2>
struct GeometryTreeWalker : public ReductionResult {
   
    typedef typename Kernel::Scalar Scalar;
    
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Parameter; //the parameters needed to describe the current dof
    G1<Kernel, false> geometry1;
    G2<Kernel, false> geometry2;
    

    Node<Kernel, G1, G2>* ResultNode; //the node we stopped at
    //std::vector<> ConstraintPool; //all remaining constraints
};

template<typename Final>
struct EdgeReductionTree {
    
    virtual void apply(typename Final::Graph& g, graph::LocalEdge e) = 0;
};

template<typename Final, template<class, bool> class G1, template<class, bool> class G2>
struct GeometryEdgeReductionTree : public EdgeReductionTree<Final>, public GeometryNode<typename Final::Kernel, G1, G2> {
    
    typedef typename Final::Kernel  Kernel;
    typedef typename Kernel::Scalar Scalar;
    
    virtual void apply(typename Final::Graph& g, graph::LocalEdge e) {
        
        typedef typename Final::Graph::global_edge_iterator iterator;
        std::pair<iterator, iterator> it = g.getGlobalEdges(e);
        
        //extract the geometry data
        symbolic::Geometry* g1 = g.template getProperty<symbolic::GeometryProperty>(boost::source(e, g));
        symbolic::Geometry* g2 = g.template getProperty<symbolic::GeometryProperty>(boost::target(e, g));
        
        dcm_assert(g1 != nullptr && g2 != nullptr);
        
        //order the symbolic constraints to match the template arguments
        symbolic::Geometry* sg1 = (Final::template geometryIndex<G1>::value == g1->type) ? g1 : g2;
        symbolic::Geometry* sg2 = (Final::template geometryIndex<G1>::value == g2->type) ? g2 : g1;
        
        dcm_assert(sg1 != sg2); 
        
        //get the primitive geometries
        const G1<Kernel, false>& tg1 = static_cast<TypeGeometry<Kernel, G1>*>(sg1)->getPrimitveGeometry();
        const G2<Kernel, false>& tg2 = static_cast<TypeGeometry<Kernel, G2>*>(sg2)->getPrimitveGeometry();
        
        //get or create and apply a treewalker which holds the data and the results
        ReductionResult* res = g.template getProperty<symbolic::ResultProperty>(e);
        GeometryTreeWalker<Kernel, G1, G2>* walker;
        if(!res)
            walker = new GeometryTreeWalker<Kernel, G1, G2>;
        else 
            walker = static_cast<GeometryTreeWalker<Kernel, G1, G2>*>(res);
        
        walker->geometry1 = tg1;
        walker->geometry2 = tg2;
        walker->Parameter = Eigen::Matrix<Scalar, 1, 1>::Zero();
        
        //calculate and store the result
        GeometryNode<Kernel, G1, G2>::applyWalker(walker);        
        g.template setProperty<symbolic::ResultProperty>(e, walker);
    }; 
    
    void execute(symbolic::GeometryTreeWalker<Kernel, G1, G2>* walker) const {
        std::cout<<"excute tree"<<std::endl;
    };
};

}//symbolic    
}//dcm

#endif //DCM_ANALYSE_H

