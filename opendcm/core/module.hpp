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
#include "reduction.hpp"
#include "constraint.hpp"
#include "geometry.hpp"
#include "solver.hpp"

#include <boost/preprocessor/seq/for_each_i.hpp>
#include <boost/preprocessor/seq/enum.hpp>
#include <boost/preprocessor/seq/size.hpp>
#include <boost/preprocessor/seq/fold_left.hpp>
#include <boost/preprocessor/seq/transform.hpp>

#include <boost/mpl/int.hpp>

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
  
/**************************************/

#define QUALIFY(s, data, elem) \
    elem<typename data::Kernel>

#define DCM_MODULE_ADD_GEOMETRIES(stacked, seq) \
    typedef mpl::vector<BOOST_PP_SEQ_ENUM(BOOST_PP_SEQ_TRANSFORM(QUALIFY, stacked, seq))> TmpGeometryList;\
    typedef typename mpl::fold<TmpGeometryList, typename stacked::GeometryList, \
         mpl::push_back<mpl::_1, mpl::_2>>::type GeometryList;
   
#define DCM_MODULE_ADD_CONSTRAINTS(stacked, seq) \
    typedef mpl::vector<BOOST_PP_SEQ_ENUM(seq)> TmpConstraintList;\
    typedef typename mpl::fold<TmpConstraintList, typename stacked::ConstraintList, \
        mpl::push_back<mpl::_1, mpl::_2>>::type ConstraintList;

namespace dcm {
namespace details {

template<typename List, typename Obj>
struct remove_base {
    typedef typename mpl::remove_if<List, mpl::and_<boost::is_base_of<mpl::_1, Obj>,
            mpl::not_<boost::is_same<mpl::_1, Obj> > > >::type type;
};


template<typename Final, typename MathKernel>
struct ModuleCoreInit {

    ModuleCoreInit() : m_graph(NULL) {
#ifdef DCM_USE_LOGGING
        sink = init_log();
        log.add_attribute("Tag",  attrs::constant<std::string>("System"));
#endif
    };

    ~ModuleCoreInit() {
#ifdef DCM_USE_LOGGING
        stop_log(sink);
#endif
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
    typedef details::Object<Final>      Object;
    typedef mpl::vector0<>              FullObjectList;
    
    //initialize the geometry and constraint handling
    typedef mpl::vector<>                               GeometryList;
    typedef mpl::vector3<Distance, Orientation, Angle>  ConstraintList;

protected:
    typedef mpl::vector<numeric::EquationBuilderProperty<Kernel>>           EdgeProperties;
    typedef mpl::vector<GraphObjectProperty, symbolic::ConstraintProperty>  GlobalEdgeProperties;
    typedef mpl::vector<GraphObjectProperty, symbolic::GeometryProperty>    VertexProperties;
    typedef mpl::vector0<>                                                  ClusterProperties;
   
#ifdef DCM_TESTING
public:
#endif
    //ensure that the correct graph type is used by not allowing anyone to set the graph pointer
    std::shared_ptr<graph::AccessGraphBase> getGraph() {
        if(!m_graph)
            m_graph = std::make_shared<typename Final::Graph>();

        return m_graph;
    };
    
protected:
#ifdef DCM_USE_LOGGING
    details::dcm_logger log;
    boost::shared_ptr< sink_t > sink;
#endif
private:
    std::shared_ptr<graph::AccessGraphBase> m_graph;    
};

template<typename Final, typename Stacked>
struct ModuleCoreFinish : public Stacked {
    
    //the math kernel in use
    typedef typename Stacked::Kernel Kernel;
    
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
    struct initGeometryIndex : utilities::index<typename Stacked::GeometryList, G>{};
    
    template<template<class> class G>
    struct geometryIndex {
        typedef typename initGeometryIndex<G<Kernel>>::type type; 
        const static long value = type::value;
    };
   
    template<typename G>
    struct constraintIndex : utilities::index<typename Stacked::ConstraintList, G>{};
    
    typedef graph::ClusterGraph<typename Stacked::EdgeProperties, typename Stacked::GlobalEdgeProperties,
            typename Stacked::VertexProperties, typename Stacked::ClusterProperties> Graph;        

            
    /**
     * @brief Solves the constraint geometry system 
     * 
     * Follow the solve procedure as outlined in the uml files
     */
    void solve() {       
        
        //build up the system and solve
        auto g = std::static_pointer_cast<Graph>(this->getGraph());
#ifdef DCM_USE_LOGGING
        BOOST_LOG_SEV(Stacked::log, details::solving) << "Setup Solver";
#endif
        auto solvable = solver::Builder::createSolvableSystem<Final>(g, m_converter);
#ifdef DCM_USE_LOGGING
        BOOST_LOG_SEV(Stacked::log, details::solving) << "Execute Solver";
#endif
        solvable->execute();
                       
        //post process the finished calculation. As solving throws an exception when not successfull we 
        //are sure that we finished and that we can process the solution
                
        
    };
    
private:
    symbolic::NumericConverter<Kernel, typename Stacked::GeometryList, 
                               typename Stacked::ConstraintList, Graph> m_converter;
};


} //details
} //dcm
#endif //DCM_MODULE_CORE_H
