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
#include "constraint.hpp"
#include "scheduler.hpp"

namespace dcm {
namespace symbolic {

/**
 * @brief Structure to setup the numeric solver from reduction results
 *
 * This struct is the entry point for the numeric solver setup. It has all the functions needed to get the
 * numeric geometries and constraints.
 */
template<typename Kernel>
struct ReductionResult {
    /*
      typedef shedule::FlowGraph::Node  Node;

      virtual numeric::GeomertyEquation<Kernel>*  createGenericGeometry(symbolic::Geometry& geom) = 0;
      virtual std::pair<numeric::GeomertyEquation<Kernel>*, Node>
                                          createGeometry(symbolic::Geometry& sknown,
                                                         numeric::GeomertyEquation<Kernel>* nknown) = 0;

      virtual std::pair<std::vector<numeric::Equation<Kernel>*>, Node>
                                          setupGenericEquations(numeric::GeomertyEquation<Kernel>* g1,
                                                                numeric::GeomertyEquation<Kernel>* g2,
                                                                shedule::FlowGraph& flow) = 0;
      virtual std::pair<std::vector<numeric::Equation<Kernel>*>, Node>
                                          setupEquations(numeric::GeomertyEquation<Kernel>* g1,
                                                         numeric::GeomertyEquation<Kernel>* g2,
                                                         shedule::FlowGraph& flow) = 0;

                                                         */
};

template<typename Kernel>
struct ResultProperty {
    typedef ReductionResult<Kernel>* type;
    struct default_value {
        ReductionResult<Kernel>* operator()() {
            return nullptr;
        };
    };
};

template<typename Kernel, template<class, bool> class G1, template<class, bool> class G2>
struct TypedReductionResult : public ReductionResult<Kernel> {
    /*
    typedef typename Final::Kernel      Kernel;
    typedef shedule::FlowGraph::Node    Node;


    virtual numeric::GeomertyEquation<Kernel>* createGenericGeometry(Geometry& geom) {

        if(geom.type == Final::template primitiveGeometryIndex<G1>::value)
            return new numeric::Geometry<Kernel, G1>();

        else
            return new numeric::Geometry<Kernel, G2>();
    };

    virtual std::pair<std::vector<numeric::Equation<Kernel>*>, Node>
                                        setupGenericEquations(numeric::GeomertyEquation<Kernel>* g1,
                                                              numeric::GeomertyEquation<Kernel>* g2,
                                                              shedule::FlowGraph& flow) {

        std::vector<numeric::Equation<Kernel>*> vec;
    };
    */
protected:
    std::vector<symbolic::Constraint*>  ConstraintPool; //all available constraints
};

template<typename Kernel, template<class, bool> class G1, template<class, bool> class G2>
struct GeometryNode;

template<typename Kernel, template<class, bool> class G1, template<class, bool> class G2>
struct GeometryTreeWalker;

/**
 * @brief Connects two nodes and is used to symbolicly solve constraints
 *
 * This class describes the connection between two GeometryNodes. This means it is responsible for the
 * transformation of one dof state to annother. For this it access the user constraints and checks if a
 * dof reduction is possible. Furthermore it is responsible for a redundency and conflict check. To use this
 * class derive you custom GeometryEdge and override the \ref apply method. Return true if the processing of the
 * GeometryTreeWalker was successfull and if we shall proceed to the next node. If not return false.
 */
template<typename Kernel, template<class, bool> class G1, template<class, bool> class G2>
struct GeometryEdge {

    GeometryNode<Kernel, G1, G2>* start;
    GeometryNode<Kernel, G1, G2>* end;

    virtual bool apply(symbolic::GeometryTreeWalker<Kernel, G1, G2>* walker) const = 0;
};

/**
 * @brief Describes the dof between two geometries
 *
 * A reduction tree node describes how much and which degrees of freedom exist between two geometries. This
 * information can be used to create the appropriate numeric depending geometries.
 */
template<typename Kernel, template<class, bool> class G1, template<class, bool> class G2>
struct GeometryNode  {

    /**
     * @brief Connect this node to annother via a specific GeometryEdge
     *
     * This function can be used to connect a node to annother one. As the connections are always conditional
     * the GeometryEdge evaluating the condition needs to be supplied.
     *
     * \param node The node we build a connection to
     * \param edge The GeometryEdge which evaluates the condition for the transition
     */
    void connect(GeometryNode<Kernel, G1, G2>* node, GeometryEdge<Kernel, G1, G2>* edge) {

        edge->start = this;
        edge->end   = node;
        edges.push_back(edge);
    };

    /**
     * @brief ${...}
     *
     * @return void
     */
    void execute(symbolic::GeometryTreeWalker<Kernel, G1, G2>* walker) {

    };

    /**
     * @brief Executed to setup the numerical solver
     *
     * This function is responsible for setting up the numerical system by creating the appropriate calculation
     * nodes. The default implementation of this function creates an equation for every remaining constraint.
     * Derived classes should always call this base version to ensure that remaining constraints are handled
     * properly.
     * @param walker The tree walker local storage
     */
    virtual void buildNumeric(symbolic::GeometryTreeWalker<Kernel, G1, G2>* walker) {
        /*
        for(symbolic::Constraint* c : walker->ConstraintPool) {


        } */
    };

protected:
    bool applyWalker(GeometryTreeWalker<Kernel, G1, G2>* walker) {
        for (const GeometryEdge<Kernel, G1, G2>& e : edges) {
            if (e.apply(walker))
                break;
        };
        return true;
    };

    std::vector<GeometryEdge<Kernel, G1, G2>> edges;
};


template<typename Final, template<class, bool> class G1, template<class, bool> class G2>
struct GeometryTreeWalker : public TypedReductionResult<Final, G1, G2> {
/*
    typedef typename Final::Kernel              Kernel;
    typedef typename Kernel::Scalar             Scalar;

    GeometryTreeWalker(const G1<Kernel, false>& g1, const G2<Kernel, false>& g2) :
        geometry1(g1), geometry2(g2) {};

//protected:
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Parameter; //the parameters needed to describe the current dof
    const G1<Kernel, false>& geometry1;
    const G2<Kernel, false>& geometry2;

    GeometryNode<Kernel, G1, G2>*       ResultNode;     //the node we stopped at*/
};

template<typename Final>
struct EdgeReductionTree {

    virtual ~EdgeReductionTree() {};

    /**
     * @brief Analyses the global edges and finds the best reduction result
     *
     * @remark The function is reentrant but is not safe to be called on the same data from multiple threads
     *
     * @param g The graph the local edge belongs to
     * @param e The local edge to analyse for reduction
     * @return void
     */
    virtual void apply(std::shared_ptr<typename Final::Graph> g, graph::LocalEdge e) = 0;
};

template<typename Final, template<class, bool> class G1, template<class, bool> class G2>
struct GeometryEdgeReductionTree : public EdgeReductionTree<Final>, public GeometryNode<typename Final::Kernel, G1, G2> {

    typedef typename Final::Kernel  Kernel;
    typedef typename Kernel::Scalar Scalar;

    virtual ~GeometryEdgeReductionTree() {};

    template<typename T>
    void pretty(const T& t) {
        std::cout<<__PRETTY_FUNCTION__<<std::endl;
    };

    virtual void apply(std::shared_ptr<typename Final::Graph> g, graph::LocalEdge e) {

/*
        //extract the geometry data
        symbolic::Geometry* g1 = g->template getProperty<symbolic::GeometryProperty>(g->source(e));
        symbolic::Geometry* g2 = g->template getProperty<symbolic::GeometryProperty>(g->target(e));

        //order the symbolic constraints to match the template arguments
        symbolic::Geometry* sg1 = (Final::template primitiveGeometryIndex<G1>::value == g1->type) ? g1 : g2;
        symbolic::Geometry* sg2 = (Final::template primitiveGeometryIndex<G1>::value == g2->type) ? g2 : g1;

        dcm_assert(sg1 != sg2);

        //get the primitive geometries
        const G1<Kernel, false>& tg1 = static_cast<TypeGeometry<Kernel, G1>*>(sg1)->getPrimitveGeometry();
        const G2<Kernel, false>& tg2 = static_cast<TypeGeometry<Kernel, G2>*>(sg2)->getPrimitveGeometry();

        //get or create a treewalker which holds the data and the results
        ReductionResult<Kernel>* res = g->template getProperty<symbolic::ResultProperty<Kernel>>(e);
        GeometryTreeWalker<Kernel, G1, G2>* walker;
        if (!res)
            walker = new GeometryTreeWalker<Kernel, G1, G2>(tg1, tg2);
        else {
            walker = static_cast<GeometryTreeWalker<Kernel, G1, G2>*>(res);
            walker->ConstraintPool.clear();
        }

        walker->ResultNode = this; //ensure default constraint handling when no reduction is possible
        walker->Parameter  = Eigen::Matrix<Scalar, 1, 1>::Zero(); //setup default values

        //setup the ConstraintPool
        typedef typename Final::Graph::global_edge_iterator iterator;
        std::pair<iterator, iterator> it = g->getGlobalEdges(e);
        for (; it.first != it.second; ++it.first)
            walker->ConstraintPool.push_back(g->template getProperty<symbolic::ConstraintProperty>(*it.first));

        //calculate and store the result
        GeometryNode<Kernel, G1, G2>::applyWalker(walker);
        g->template setProperty<symbolic::ResultProperty>(e, walker);*/
    };
};

}//symbolic
}//dcm

#endif //DCM_ANALYSE_H

