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

#ifndef DCM_MODULE_CORE_H
#define DCM_MODULE_CORE_H

#include "accessgraph.hpp"
#include "clustergraph.hpp"
#include "filtergraph.hpp"
#include "logging.hpp"
#include "object.hpp"
#include "analyse.hpp"
#include "constraint.hpp"
#include "geometry.hpp"
#include "solver.hpp"

#include <boost/preprocessor/seq/for_each_i.hpp>
#include <boost/preprocessor/seq/enum.hpp>
#include <boost/preprocessor/seq/size.hpp>
#include <boost/preprocessor/seq/fold_left.hpp>
#include <boost/preprocessor/seq/transform.hpp>

#include <boost/mpl/int.hpp>
#include <boost/multi_array.hpp>


namespace mpl = boost::mpl;

#define DCM_MODULE_ADD_OBJECTS(stacked, seq) \
    typedef mpl::vector<BOOST_PP_SEQ_ENUM(seq)> TmpFullObjectList;\
    typedef typename mpl::fold<TmpFullObjectList, typename stacked::FullObjectList, \
        mpl::push_back<mpl::_1, mpl::_2>>::type FullObjectList;

#define DCM_MODULE_ADD_VERTEX_PROPERTIES(stacked, seq) \
    typedef mpl::vector<BOOST_PP_SEQ_ENUM(seq)> TmpVertexProperties;\
    typedef typename mpl::fold<TmpVertexProperties, typename stacked::VertexProperties, \
         mpl::push_back<mpl::_1, mpl::_2>>::type VertexProperties;
         
#define DCM_MODULE_ADD_LOCAL_EDGE_PROPERTIES(stacked, seq) \
    typedef mpl::vector<BOOST_PP_SEQ_ENUM(seq)> TmpEdgeProperties;\
    typedef typename mpl::fold<TmpEdgeProperties, typename stacked::EdgeProperties, \
         mpl::push_back<mpl::_1, mpl::_2>>::type EdgeProperties;

#define DCM_MODULE_ADD_GLOBAL_EDGE_PROPERTIES(stacked, seq) \
    typedef mpl::vector<BOOST_PP_SEQ_ENUM(seq)> TmpGlobalEdgeProperties;\
    typedef typename mpl::fold<TmpGlobalEdgeProperties, typename stacked::GlobalEdgeProperties, \
         mpl::push_back<mpl::_1, mpl::_2>>::type GlobalEdgeProperties;
         
#define DCM_MODULE_ADD_CLUSTER_PROPERTIES(stacked, seq) \
    typedef mpl::vector<BOOST_PP_SEQ_ENUM(seq)> TmpClusterProperties;\
    typedef typename mpl::fold<TmpClusterProperties, typename stacked::ClusterProperties, \
         mpl::push_back<mpl::_1, mpl::_2>>::type ClusterProperties;
         

namespace dcm {
namespace details {

template<typename List, typename Obj>
struct remove_base {
    typedef typename mpl::remove_if<List, mpl::and_<boost::is_base_of<mpl::_1, Obj>,
            mpl::not_<boost::is_same<mpl::_1, Obj> > > >::type type;
};


template<typename Final, typename MathKernel>
struct ModuleCoreInit {

    ModuleCoreInit() : graph(NULL)
#ifdef DCM_USE_LOGGING
        , sink(init_log())
#endif
    {
        int size = mpl::size<typename Final::GeometryList>::type::value;
        reduction.resize(boost::extents[size][size]);
        
        int constraints = mpl::size<typename Final::ConstraintList>::type::value;
        generator.resize(boost::extents[size][size][constraints]);
        
        //build up the default reduction nodes
        //mpl trickery to get a sequence counting from 0 to the size of stroage entries
        typedef mpl::range_c<int,0,
                mpl::size<typename Final::GeometryList>::value> StorageRange;
        //now iterate that sequence so we can access all storage elements with knowing the position
        //we are at
        RecursiveSequenceApplyer<typename Final::GeometryList, ReductionNodeCreator> r(reduction);
        mpl::for_each<StorageRange>(r);
        
        RecursiveSequenceApplyer<typename Final::GeometryList, ConstraintGeneratorCreator> g(generator);
        mpl::for_each<StorageRange>(g);
        
    };

    ~ModuleCoreInit() {
#ifdef DCM_USE_LOGGING
        stop_log(sink);
#endif
        /*
        //delete all reduction nodes
        int size = std::pow(mpl::size<typename Final::GeometryList>::value, 2);
        for(int i=0; i<size; ++i)
            delete reduction(i);*/
    };

#ifdef DCM_USE_LOGGING
    template<typename Expr>
    void setLoggingFilter(const Expr& ex) {
        sink->set_filter(ex);
    }
#endif

    //the math Kernel for everyone accessible and the settings hanling needed for it
    typedef MathKernel Kernel;

    //initialize the object handling
    typedef Object<Final>  ObjectBase;
    typedef mpl::vector0<> FullObjectList;
    
    //initialize the geometry and constraint handling
    typedef mpl::vector<>                               GeometryList;
    typedef mpl::vector3<Distance, Orientation, Angle>  ConstraintList;

protected:
    typedef mpl::vector<symbolic::ResultProperty<Kernel>>       EdgeProperties;
    typedef mpl::vector<symbolic::ConstraintProperty>           GlobalEdgeProperties;
    typedef mpl::vector<symbolic::GeometryProperty>             VertexProperties;
    typedef mpl::vector0<>                                      ClusterProperties;

    template<template<class, bool> class G1, template<class, bool> class G2>
    symbolic::GeometryNode<Kernel, G1, G2>* getInitialReductionNode() {
        int n1 = Final::template primitiveGeometryIndex<G1>::value;
        int n2 = Final::template primitiveGeometryIndex<G2>::value;
        return reinterpret_cast<symbolic::GeometryNode<Kernel, G1, G2>*>(reduction[n1][n2]);
    };
    
    template<template<class, bool> class G1, template<class, bool> class G2>
    symbolic::EdgeReductionTree<Final>* getInitialReductionTree() {
        int n1 = Final::template primitiveGeometryIndex<G1>::value;
        int n2 = Final::template primitiveGeometryIndex<G2>::value;
        return reduction[n1][n2];
    };
    
    template<template<class, bool> class G1, template<class, bool> class G2, typename PC>
    numeric::ConstraintEquationGenerator<Kernel> getEquationGenerator() {
        int n1 = Final::template primitiveGeometryIndex<G1>::value;
        int n2 = Final::template primitiveGeometryIndex<G2>::value;
        int n3 = Final::template constraintIndex<PC>::value;
        return reduction[n1][n2][n3];
    };
    
#ifdef DCM_TESTING
public:
#endif
    //ensure that the correct graph type is used by not allowing anyone to set the graph pointer
    std::shared_ptr<graph::AccessGraphBase> getGraph() {
        if(!graph)
            graph = std::make_shared<typename Final::Graph>();

        return graph;
    };
    
protected:
    boost::multi_array<symbolic::EdgeReductionTree<Final>*,2>           reduction;
    boost::multi_array<numeric::ConstraintEquationGenerator<Kernel>,3>  generator;

private:
    std::shared_ptr<graph::AccessGraphBase> graph;
#ifdef DCM_USE_LOGGING
    boost::shared_ptr< sink_t > sink;
#endif
    
    template<typename Sequence, template<class> class Functor>
    struct RecursiveSequenceApplyer {
        
        Functor<Sequence> functor;
        
        template<typename T>
        RecursiveSequenceApplyer(T& param) 
            : functor(Functor<Sequence>(param)) {};
            
        RecursiveSequenceApplyer<RecursiveSequenceApplyer<Sequence, Functor>>(RecursiveSequenceApplyer& r) 
            : functor(r.functor) {};
            
        template<typename T>
        void operator()(const T& t) {        
              typedef mpl::range_c<int, T::value, mpl::size<Sequence>::value> StorageRange;
              mpl::for_each<StorageRange>(InnerLoop<T>(functor));
        };
        
        template<typename Number>
        struct InnerLoop {
    
            Functor<Sequence>& functor;
            
            InnerLoop(Functor<Sequence>& f) 
                : functor(f) {};
                
            template<typename T>
            void operator()(const T& t) {            
                functor.template operator()<Number, T>();
            };
        };
    };
    
    template<typename Sequence>
    struct ReductionNodeCreator {
    
        boost::multi_array<symbolic::EdgeReductionTree<Final>*,2>& reduction;
        
        ReductionNodeCreator(boost::multi_array<symbolic::EdgeReductionTree<Final>*,2>& r) 
            : reduction(r) {};
            
        template<typename N1, typename N2>
        void operator()() {
        
            typedef typename mpl::at<Sequence, N1>::type t1;
            typedef typename mpl::at<Sequence, N2>::type t2;
            
            auto node = new symbolic::GeometryEdgeReductionTree<Final, 
                                geometry::extractor<t1>::template primitive,
                                geometry::extractor<t2>::template primitive >();
        
            int idx1 = Final::template geometryIndex<t1>::value;
            int idx2 = Final::template geometryIndex<t2>::value;
            
            reduction[idx1][idx2] = node;
            reduction[idx2][idx1] = node;
        };
    };
    
    template<typename Sequence>
    struct ConstraintGeneratorCreator {
    
        boost::multi_array<numeric::ConstraintEquationGenerator<Kernel>,3>& generator;
        
        ConstraintGeneratorCreator(boost::multi_array<numeric::ConstraintEquationGenerator<Kernel>,3>& r) 
            : generator(r) {};
            
        template<typename N1, typename N2>
        void operator()() {
        
            
        };
        
        template<typename G1, typename G2, int n1, int n2>
        struct InnerLoop {
            
            boost::multi_array<numeric::ConstraintEquationGenerator<Kernel>,3>& generator;
        
            InnerLoop(boost::multi_array<numeric::ConstraintEquationGenerator<Kernel>,3>& r) : generator(r) {};
            
            template<typename T>
            void operator()(const T& t) {
            
                generator[n1][n2][T::value] = numeric::TypedConstraintEquationGenerator<Kernel, 
                                                  typename mpl::at<typename Final::ConstraintList, T>::type,
                                                  geometry::extractor<G1>::template primitive,
                                                  geometry::extractor<G2>::template primitive>();
                                                                         
                generator[n2][n1][T::value] = generator[n1][n2][T::value];
            };
        };
    };
};

template<typename Final, typename Stacked>
struct ModuleCoreFinish : public Stacked {
    
    //The FullObjectList holds all created object types, including all "evolution steps", meaning every
    //base class of the final object is represented as well. However, we often want only to habe the user-visible types
    //without the base classes. Therefore the ObjectList is provided which only holds types which are directly
    //visible to the user.
    typedef typename mpl::fold<typename Stacked::FullObjectList, typename Stacked::FullObjectList,
            remove_base< mpl::_1, mpl::_2 > >::type ObjectList;

    //Metafunction to filter out the exact object type which holds a specific property. As a property can by definition
    //only added once to a object this operation should yield the only relevant type.
    template<typename Prop>
    struct objectByProperty {
        typedef typename mpl::find_if<typename Stacked::FullObjectList, object_has_property<mpl::_1, Prop> >::type iter;
        BOOST_MPL_ASSERT((mpl::not_< boost::is_same<iter, typename mpl::end<ObjectList>::type> >));

        typedef typename mpl::deref< iter >::type type;
    };

    //Metafunction to extract the objects type id as integral constant
    template<typename Object>
    struct objectTypeID {
        typedef typename mpl::find<ObjectList, Object>::type iterator;
        typedef typename mpl::distance<typename mpl::begin<ObjectList>::type, iterator>::type ID;
    };
    
    template<typename G>
    struct geometryIndex {
        typedef typename mpl::find<typename Stacked::GeometryList, G>::type iterator;
        typedef boost::is_same<iterator, typename mpl::end<typename Stacked::GeometryList>::type> valid;
        BOOST_MPL_ASSERT_MSG(mpl::not_<valid>::value, GEOMETRY_TYPE_NOT_REGISTERT, (G));
        
        typedef typename mpl::distance<typename mpl::begin<typename Stacked::GeometryList>::type,
                            iterator>::type type; 
        const static long value = type::value;
    };

    template<template<class, bool> class G>
    struct geometryIndex<geometry::adaptor<G>> {
        typedef typename geometryIndex<typename geometry::adaptor<G>::placeholder>::type type; 
        const static long value = type::value;
    };
    
    template<template<class, bool> class G>
    struct primitiveGeometryIndex {
        typedef typename geometryIndex<typename geometry::adaptor<G>::placeholder>::type type; 
        const static long value = type::value;
    };
    
    template<typename G>
    struct constraintIndex {
        typedef typename mpl::find<typename Stacked::ConstraintList, G>::type iterator;
        typedef boost::is_same<iterator, typename mpl::end<typename Stacked::ConstraintList>::type> valid;
        BOOST_MPL_ASSERT_MSG(mpl::not_<valid>::value, CONSTRAINT_TYPE_NOT_REGISTERT, (G));
        
        typedef typename mpl::distance<typename mpl::begin<typename Stacked::ConstraintList>::type,
                            iterator>::type type; 
        const static long value = type::value;
    };
    
    typedef graph::ClusterGraph<typename Stacked::EdgeProperties, typename Stacked::GlobalEdgeProperties,
            typename Stacked::VertexProperties, typename Stacked::ClusterProperties> Graph;
         
    using Stacked::reduction;
    
        

            
    /**
     * @brief Solves the constraint geometry system 
     * 
     * Follow the solve procedure as outlined in the uml files
     */
    void solve() {       
        
        //All graph manipulation work has been done, from here on we only access the graph. 
        //Next find all connected components and build the numeric solving system based on them
        solver::createSolvableSystem(std::static_pointer_cast<Graph>(this->getGraph()), reduction);
                
        //post process the finished calculation
        
    };
};


} //details
} //dcm
#endif //DCM_MODULE_CORE_H
