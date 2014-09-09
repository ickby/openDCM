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

#include <boost/mpl/int.hpp>

#include "clustergraph.hpp"
#include "logging.hpp"
#include "object.hpp"

#include <boost/preprocessor/seq/for_each_i.hpp>
#include <boost/preprocessor/seq/enum.hpp>
#include <boost/preprocessor/seq/size.hpp>
#include <boost/preprocessor/seq/fold_left.hpp>

#define DCM_MODULE_ADD_OBJECTS(stacked, seq) \
    typedef mpl::vector<BOOST_PP_SEQ_ENUM(seq)> TmpFullObjectList;\
    typedef typename mpl::fold<TmpFullObjectList, typename stacked::FullObjectList, \
        mpl::push_back<mpl::_1, mpl::_2> >::type FullObjectList;

namespace dcm {
namespace details {

namespace mpl = boost::mpl;

template<typename List, typename Obj>
struct remove_base {
    typedef typename mpl::remove_if<List, mpl::and_<boost::is_base_of<mpl::_1, Obj>,
            mpl::not_<boost::is_same<mpl::_1, Obj> > > >::type type;
};


template<typename Final>
struct ModuleCoreInit {

    ModuleCoreInit() : graph(NULL)
#ifdef USE_LOGGING
        , sink(init_log())
#endif
    {};

    ~ModuleCoreInit() {
        if(graph)
            delete graph;

#ifdef USE_LOGGING
        stop_log(sink);
#endif
    };

#ifdef USE_LOGGING
    template<typename Expr>
    void setLoggingFilter(const Expr& ex) {
        sink->set_filter(ex);
    }
#endif

    //initialize the object handling
    typedef Object<Final>  ObjectBase;
    typedef mpl::vector0<> FullObjectList;

protected:
    typedef mpl::vector0<> EdgeProperties;
    typedef mpl::vector0<> GlobalEdgeProperties;
    typedef mpl::vector0<> VertexProperties;
    typedef mpl::vector0<> ClusterProperties;

    //ensure that the correct graph type is used by not allowing anyone to set the graph pointer
    ClusterGraphBase* getGraph() {
        if(!graph)
            graph = new typename Final::Graph();

        return graph;
    };

private:
    ClusterGraphBase* graph;
#ifdef USE_LOGGING
    boost::shared_ptr< sink_t > sink;
#endif

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

protected:
    typedef ClusterGraph<typename Stacked::EdgeProperties, typename Stacked::VertexProperties,
            typename Stacked::ClusterProperties, typename Stacked::GlobalEdgeProperties> Graph;
};

} //details
} //dcm
#endif //DCM_MODULE_CORE_H
