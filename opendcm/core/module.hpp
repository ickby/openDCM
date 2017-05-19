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

#include <chrono>

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
   
#define DCM_MODULE_ADD_BINARY CONSTRAINTS(stacked, seq) \
    typedef mpl::vector<BOOST_PP_SEQ_ENUM(seq)> TmpBinaryConstraintList;\
    typedef typename mpl::fold<TmpBinaryConstraintList, typename stacked::BinaryConstraintList, \
        mpl::push_back<mpl::_1, mpl::_2>>::type BinaryConstraintList;

#define DCM_MODULE_ADD_UNARY_CONSTRAINTS(stacked, seq) \
    typedef mpl::vector<BOOST_PP_SEQ_ENUM(seq)> TmpUnaryConstraintList;\
    typedef typename mpl::fold<TmpUnaryConstraintList, typename stacked::UnaryConstraintList, \
        mpl::push_back<mpl::_1, mpl::_2>>::type UnaryConstraintList;
        
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
        
        //ensure the graph and fix cluster is initialised
        getFixCluster();
    };

    ~ModuleCoreInit() {
#ifdef DCM_USE_LOGGING
        stop_log(sink);
#endif
    };
    
    //this function allows to initialize modules after all constructors are done. It wil be called 
    //for every module that implements it, and it is than that modules responsibility to call init 
    //on its base class.
    void init() {};

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
    typedef mpl::vector3<Distance, Orientation, Angle>  BinaryConstraintList;
    typedef mpl::vector1<Fix>                           UnaryConstraintList;
    
protected:
    typedef mpl::vector1<numeric::EquationHandlerProperty<Kernel>> EdgeProperties;
    typedef mpl::vector2<GraphObjectProperty, 
                         symbolic::ConstraintProperty>             GlobalEdgeProperties;
    typedef mpl::vector3<GraphObjectProperty, 
                         symbolic::GeometryProperty, 
                         symbolic::ConstraintListProperty>         VertexProperties;
    typedef mpl::vector0<>                                         ClusterProperties;
   
#ifdef DCM_USE_LOGGING
    details::dcm_logger log;
    boost::shared_ptr< sink_t > sink;
#endif
    
    
    //function for setting up the reduction tree in the modules. As we don't know the generator 
    //type of tree types now we simply go with the source node, that is all that is needed
    virtual reduction::EdgeReductionGraph* getReductionGraph(int Type1, int Type2) = 0;

#ifdef DCM_TESTING
public:
#endif
    //ensure that the correct graph type is used by not allowing anyone to set the graph pointer
    std::shared_ptr<graph::AccessGraphBase> getGraph() {
        if(!m_graph)
            m_graph = std::make_shared<typename Final::Graph>();

        return m_graph;
    };
    
    std::shared_ptr<graph::AccessGraphBase> getFixCluster() {
        
        if(!m_fixCluster) {
            std::shared_ptr<typename Final::Graph> cluster = std::static_pointer_cast<typename Final::Graph>(this->getGraph());  
            auto result = cluster->createCluster();
            m_fixCluster = result.first;
            result.first->template setProperty<graph::Type>(graph::FixCluster);
            cluster->template setProperty<graph::Type>(result.second, graph::FixCluster);
        }
            
        return m_fixCluster;
    };
    
private:
    std::shared_ptr<graph::AccessGraphBase> m_graph, m_fixCluster;    
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


    //Metafunction to extract the objects type id as integral constant
    template<typename Object>
    struct objectTypeID {
        typedef typename mpl::find<ObjectList, Object>::type iterator;
        typedef typename mpl::distance<typename mpl::begin<ObjectList>::type, iterator>::type ID;
    };
    
    typedef graph::ClusterGraph<typename Stacked::EdgeProperties, typename Stacked::GlobalEdgeProperties,
            typename Stacked::VertexProperties, typename Stacked::ClusterProperties> Graph;        

            
    ModuleCoreFinish() {
    
        //now all relevant constructors are done, so let's init everything!
        //calling init should be processed by whatever module implements it and than be passed upwards.
        Stacked::init();
    };
            
    /**
     * @brief Solves the constraint geometry system 
     * 
     * Follow the solve procedure as outlined in the uml files
     */
    void solve() {       
        
        auto start = std::chrono::system_clock::now();
     
        auto g = std::static_pointer_cast<Graph>(this->getGraph());  
        
#ifdef DCM_USE_LOGGING
        BOOST_LOG_SEV(Stacked::log, details::solving) << "Start Preprocessing";
#endif
        //pre process to get the data into the system. 
        solver::Builder<Kernel>::preprocessSystem(g);
        
#ifdef DCM_USE_LOGGING
        BOOST_LOG_SEV(Stacked::log, details::solving) << "Execute Solver";
#endif
        solver::Builder<Kernel>::template solveSystem<Final>(g, m_converter);
        
#ifdef DCM_USE_LOGGING
        BOOST_LOG_SEV(Stacked::log, details::solving) << "Done Solving, start postprocessing";
#endif
        //post process the finished calculation. As solving throws an exception when not successfull we 
        //are sure that we finished and that we can process the solution
        solver::Builder<Kernel>::postprocessSystem(g);
#ifdef DCM_USE_LOGGING
        BOOST_LOG_SEV(Stacked::log, details::solving) << "Done postprocessing. Solver time: ";
#endif 
        auto end = std::chrono::system_clock::now();
        std::cout<<"solver time: "<<std::chrono::duration_cast<std::chrono::microseconds>(end-start).count()<<std::endl;
    };
    
protected:
    //function for setting up the reduction tree in the modules. As we don't know the generator 
    //type of tree types now we simply go with the source node, that is all that is needed
    virtual reduction::EdgeReductionGraph* getReductionGraph(int Type1, int Type2) override {
        return m_converter.getReductionGraph(Type1, Type2);
    };
    
private:
    symbolic::NumericConverter<Kernel, typename Stacked::GeometryList, typename Stacked::BinaryConstraintList,
                               typename Stacked::UnaryConstraintList, Graph> m_converter;
};


} //details
} //dcm
#endif //DCM_MODULE_CORE_H
