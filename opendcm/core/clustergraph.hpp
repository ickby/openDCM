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

#ifndef DCM_CLUSTERGRAPH_HPP
#define DCM_CLUSTERGRAPH_HPP

#include "accessgraph.hpp"

namespace mpl = boost::mpl;
namespace fusion = boost::fusion;

namespace dcm {
namespace graph {
    
using namespace details;

/** @addtogroup Core
 * @{
 * */

/**
 * @brief Generator for unique identifiers
 *
 * The universalID used to identify vertices and edges globaly need to be unique and therefore can't be
 * created at good will. This generator creates universalID's in a incremental manner and is intended to
 * to be shared between all graphs of a system, so that all created ID's are unique.
 **/
struct IDgen {
    universalID* counter;

    IDgen() {
        counter = new universalID(10);
    };
    IDgen(universalID id) {
        counter = new universalID(id);
    };
    ~IDgen() {
        delete counter;
    };
    /**
     * @brief Generates a new unique ID
     *
     * @return :universalID
     **/
    universalID generate() {
        return ++ (*counter);
    };
    /**
     * @brief Returns the amount if generated ID's
     *
     * As universalID's are integers the returned count is a ID and can therefore also be used as the last
     * created ID.
     *
     * @return :universalID
     **/
    universalID count() {
        return (*counter);
    };
    /**
     * @brief Set the current value for incremental creation
     *
     * ID's are created incrementaly and if a specific startingpoint is wished it can be set here by
     * supplying the last created ID or the amount of totaly created ID's
     *
     * @param id The last created ID
     * @return void
     **/
    void setCount(universalID id) {
        *counter = id;
    };
};

/**
 * @brief Pointer type to share a common ID generator @ref IDgen
 **/
typedef std::shared_ptr<IDgen> IDpointer;

template<typename T1, typename T2, typename T3, typename T4, typename T5>
using adjacency_list = boost::adjacency_list<T1,T2,T3,T4,T5>;

/**
 * @ingroup ClusterGraph
 * @brief A graph that can be stacked in a tree-like manner without loosing its connections
 *
 * This graph implements the \ref AccessGraph interface on a boost adjacency_list. It further
 * extend the api with possibility to change the graph like adding and removing edges/vertices.
 * Furthermore it adds support for subclusters as described in the overall graph documentation
 * \ref ClusterGraph
 *
 * @tparam edge_prop a mpl::vector with properties which are added to local edges
 * @tparam globaledge_prop a mpl::vector with properties which are added to global edges
 * @tparam vertex_prop a mpl::vector with properties which are added to vertices
 * @tparam cluster_prop a mpl::vector with properties which are added to all clusters
 **/
template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
class ClusterGraph : public AccessGraph<edge_prop, globaledge_prop, 
                                            vertex_prop, cluster_prop, adjacency_list> {

    typedef AccessGraph<edge_prop, globaledge_prop, 
                        vertex_prop, cluster_prop, adjacency_list> Base;
                                                
    typedef std::map<LocalVertex, std::shared_ptr<ClusterGraph> > ClusterMap;
    
public:
    typedef typename Base::Graph Graph;
    typedef typename Base::local_vertex_iterator local_vertex_iterator;
    typedef typename Base::local_out_edge_iterator local_out_edge_iterator;
    typedef typename Base::edge_single_iterator edge_single_iterator;
    typedef typename Base::sp_base sp_base;
    typedef typename Base::vertex_bundle vertex_bundle;
    typedef typename Base::edge_bundle_single edge_bundle_single;
    typedef typename Base::GEdgeProperty GEdgeProperty;    
    
    /**
     * @brief Iterator for clusters
     *
     * Allows to iterate over all subclusters.
     **/
    typedef typename ClusterMap::iterator cluster_iterator;
    /**
     * @brief Const equivalent to \ref cluster_iterator
     **/
    typedef typename ClusterMap::const_iterator const_cluster_iterator;
 
    /**
     * @brief Basic constructor
     *
     * This constructor creates a empty cluster with a new ID generator. This is to be used on initial
     * clustergraph creation, so only for the very first cluster.
     **/
    ClusterGraph() : Base(m_graph), m_id(new IDgen) {};

    /**
     * @brief Dependent constructor
     *
     * This constructor creates a new cluster, but uses the given cluster as parent. It will therefore
     * create a tree-like relationship. Be aware, that the new cluster is not added to the parents
     * subcluster list, that has to be done manualy. The new cluster shares the parents ID generator.
     *
     * @param g the parent cluster graph
     **/
    ClusterGraph(std::shared_ptr<ClusterGraph> g) : Base(m_graph), m_parent(g), m_id(g->m_id) {};
    
    ~ClusterGraph() {};  
    
    /**
     * @brief Copys the Clustergraph into a new one
     *
     * Copys this cluster and all subclusters into the give one, which is cleared bevore copying. Be
     * aware that all properties are only copied, and as some may hold shared pointers you may have to
     * clone them. If needed this can be done with the supplied functor, which receives all copied objects
     * to his function operator which returns the new value.
     * @param into The Graph that should be a copy of this
     * @param functor The function objects which gets the graph properties and returns the ones for the
     * copied graph
     */
    template<typename Functor>
    void copyInto(std::shared_ptr<ClusterGraph> into, Functor& functor) const;
       
    /**
     * @brief Set diffrent behaviour for changed markers
     *
     * Some methods of the AccessGraph set it's changed_prop to true. Thats sensible, as they change
     * the graph. However, there are situations where you want to use the methods but don't want the change
     * marked. For example recreations while cloning. This method can be used to disable the changed setting.
     * @param on Turn change markers on or of
     * @return void
     **/
    void setCopyMode(bool on);
 
    /**
     * @brief Set the property of a owned cluster
     *
     * Makes it easy to set a property of a subcluster without retrieving it first
     *
     * @tparam P the property type which shall be set
     * @param v the local vertex which describes the subcluster
     **/
    template<typename P>
    typename P::type& getSubclusterProperty(LocalVertex v);

    
    /* *******************************************************
     * Subclustering
     * *******************************************************/

    /**
     * @brief Creates a new subcluster
     *
     * As clusters can be stacked in a tree like manner, this function can be used to create new
     * children. It automaticly adds it to the subcluster list and adds it to the graph. The new
     * subcluster is fully defined by its object and the vertex descriptor which is it's position
     * in the current cluster.
     *
     * @return :pair< std::shared_ptr< ClusterGraph >, LocalVertex > Subcluster and its descriptor
     **/
    std::pair<std::shared_ptr<ClusterGraph>, LocalVertex> createCluster();

    /**
     * @brief Returns the parent cluster
     *
     * In the stacked cluster hirarchy most clusters have a parent whcih can be accessed with this function.
     * However, the toplevel cluster dos nothave a parent and a empty shared_ptr is returned.
     *
     * @return :shared_ptr< ClusterGraph > the parent cluster or empty pointer
     **/
    std::shared_ptr<ClusterGraph> parent();

    /**
     * @brief const version of \ref parent()
     *
     * @return :shared_ptr< ClusterGraph >
     **/
    const std::shared_ptr<ClusterGraph> parent() const;

    /**
     * @brief Is this the toplevel cluster?
     *
     * @return bool if it is
     **/
    bool isRoot() const;

    /**
     * @brief Returns the toplevel cluster
     *
     * @return :shared_ptr< ClusterGraph >
     **/
    std::shared_ptr<ClusterGraph>	 root();

    /**
     * @brief const equivalent of \ref root()
     *
     * @return :shared_ptr< ClusterGraph >
     **/
    const std::shared_ptr<ClusterGraph> root() const;

    /**
     * @brief Iterators for all subclusters
     *
     * A pair with two \ref cluster_iterator is returned which point to the first cluster and
     * to one after the last. #this allows full iteration over all subclusters
     *
     * @return :pair< cluster_iterator, cluster_iterator >
     **/
    std::pair<cluster_iterator, cluster_iterator> clusters();

    /**
     * @brief const equivalent to \ref clusters()
     *
     * @return :pair< const_cluster_iterator, const_cluster_iterator >
     **/
    std::pair<const_cluster_iterator, const_cluster_iterator> clusters() const;

    /**
     * @brief The amount of all subclusters
     *
     * @return :size_t
     **/
    std::size_t numClusters() const;

    /**
     * @brief Check if this vertex is a cluster
     *
     * A subcluster is added as normal vertex to the parent cluster. There is no way to distinguish
     * between clusters and normal vertices with global or local descriptors only. Therefore this
     * function can be used to get information about the type. If it is a cluster, it can be accessed
     * with \ref getVertexCluster
     *
     * @param v The vertex to be checked
     * @return bool is cluster or not
     **/
    bool isCluster(const LocalVertex v) const;

    /**
     * @brief Get the cluster corresponding the discriptor
     *
     * A subcluster is added as normal vertex to the parent cluster. There is no way to access
     * the clusters object with global or local descriptors only. Therefore this
     * function can be used to get the object belonging to the descriptor. If the vertex is not
     * a cluster an empty pointer is returned.
     *
     * @param v The vertex for which the cluster is wanted
     * @return std::shared_ptr<ClusterGraph> the coresponding cluster orempty pointer
     **/
    std::shared_ptr<ClusterGraph> getVertexCluster(LocalVertex v);

    /**
     * @brief Get the vertex descrptor which descripes the clusters position in the graph
     *
     * This function is the inverse to \ref getVertexCluster
     *
     * @param g the graph for which the vertex is searched
     * @return :LocalVertex
     **/
    LocalVertex	getClusterVertex(std::shared_ptr<ClusterGraph> g);

    /**
     * @brief Convinience function for \ref removeCluster
     **/
    template<typename Functor>
    void removeCluster(std::shared_ptr<ClusterGraph> g, Functor& f);
    /**
     * @brief Convinience function for \ref removeCluster
     **/
    void removeCluster(std::shared_ptr<ClusterGraph> g);
    /**
     * @brief Delete all subcluster
     *
     * @return void
     **/

    void clearClusters();

    /**
     * @brief Remove a subcluster and applys the functor to all removed edges and vertices
     *
     * All downstream elements of the local vertex v will be removed after the functor is applied to there
     * edges and vertices. Note that the LocalVertex which represents the cluster to delete is not passed
     * to the functor. When ever the cluster is changed it will be passed to the functor, so that it need
     * to have three overloads: operator()(GlobalEdge), operator()(GlobalVertex), operator()(ClusterGraph&)
     *
     * @param v Local vertex which is a cluster and which should be deleted
     * @param f Functor to apply on all graph elements
     */
    template<typename Functor>
    void removeCluster(LocalVertex v, Functor& f);
    void removeCluster(LocalVertex v);

protected:
    template<typename Functor>
    void remove_vertices(Functor& f, bool recursive = false);


    /* *******************************************************
    * Creation Handling
    * *******************************************************/

public:
    /**
     * @brief Add a vertex to the local cluster
     *
     * @return fusion::vector<LocalVertex, GlobalVertex> the local and global vertex descriptor
     **/
    fusion::vector<LocalVertex, GlobalVertex> addVertex();

    /**
     * @brief Add a vertex to the local cluster with given global identifier
     *
     * Sometimes it is needed to add a vertex with given global identifier. As the global vertex can not
     * be changed after creation, this method can be used to specify the global vertex by which this
     * graph vertex can be identified. The given global vertex is not checked, you need to ensure that
     * it is a unique id or the already existing vertex is returned.
     * The ID generator is changed so that it creates only identifier bigger than v.
     *
     * @return fusion::vector<LocalVertex, GlobalVertex> the local and global vertex descriptor
     **/
    fusion::vector<LocalVertex, GlobalVertex> addVertex(GlobalVertex v);


    /**
     * @brief Add a edge between two vertices, defined by local descriptors.
     *
     * Add an edge that connects the two vertices and in the local clustergraph and assign the GlobalEdge to it. The
     * LocalVertex parameters should not represent a cluster which would result in the functions failure. If there's
     * already a local edge between the vertices a new global edge will be added and returned. Failure will be
     * recocnisable by a false value in the returned type sequence.
     *
     * @param source The first vertex the edge should connect
     * @param target The second vertex the edge should connect
     * @return fusion::vector<LocalEdge, GlobalEdge, success> with the local and global descriptors of the edge and an bool
     * value indicationg the successful creation.
     **/
    fusion::vector<LocalEdge, GlobalEdge, bool> addEdge(LocalVertex source, LocalVertex target);

    /**
     * @brief Add a edge between two vertices, defined by global descriptors.
     *
     * Adds an edge between vertices which are not nesseccarily in this local cluster and have therefore to be
     * identified with global descriptors. The only condition for source and target vertex is that both must be
     * in the local cluster or any of its subclusters. If thats not the case, the function will fail. On success
     * a new GlobalEdge will be created, but not neccessarily a local one. If the vertices are in different cluster
     * which are already connected the global edge will be added to this connecting local edge. Thats the one returned
     * in the seqence. Note that it's possible that the local edge belongs to another subcluster and therefore can't be
     * used in the local cluster. This case is indicated by the scope return value.
     *
     * @param source The first vertex the edge should connect
     * @param target The second vertex the edge should connect
     * @return fusion:vector< LocalEdge, GlobalEdge, success, scope > with the new global edge descriptor and the local
     * one where it was added. Success indicates if the function was successful and scope shows the validy of the local
     * descriptor in this cluster (true means the edge is in this cluster).
     **/
    fusion::vector<LocalEdge, GlobalEdge, bool, bool> addEdge(GlobalVertex source, GlobalVertex target);

    fusion::vector<LocalEdge, GlobalEdge, bool, bool> addEdgeGlobal(GlobalVertex source, GlobalVertex target);

    /**
     * @brief Get the local edge which holds the specified global one and the subcluster in which it is valid.
     *
     * The function only fails when the global edge is hold by a local one upstream in the cluster
     * herarchy.
     *
     * @param e GlobalEdge for which the containing local one is wanted
     * @return fusion::vector<LocalEdge, ClusterGraph*, bool> with the containing LocalEdge, the cluster which holds it and a bool indicator if function was successful.
     **/
    fusion::vector<LocalEdge, ClusterGraph*, bool> getLocalEdgeGraph(GlobalEdge e);

     /**
     * @brief Get the local vertex which holds the specified global one and the subcluster in which it is valid.
     *
     * The function only fails when the global vertex is hold by a local one upstream in the cluster
     * herarchy.
     *
     * @param v GlobalVertex for which the containing local one is wanted
     * @return fusion::vector<LocalVertex, ClusterGraph*, bool> with the containing LocalVertex, the cluster which holds it and a bool indicator if function was successful.
     **/
    fusion::vector<LocalVertex, std::shared_ptr<ClusterGraph>, bool> getLocalVertexGraph(GlobalVertex v);


    /* *******************************************************
     * Remove Handling
     * *******************************************************/
private:

    template<typename Functor>
    void downstreamRemoveVertex(GlobalVertex v, Functor& f);

    void simpleRemoveEdge(LocalEdge e);


public:
    /**
    * @brief Removes a vertex from the local cluster and applys functor to removed edges
    *
    * Removes the vertex from the local graph and invalidates the global vertex id. Also all edges connecting
    * to this vertex will be removed after the functor was applied to them. The functor needs to implement
    * operato()(GlobalEdge e). Remark that there is no checking done if the vertex is a cluster, so you
    * need to make sure it's not, as removing a clustervertex will not delete the coresponding cluster.
    *
    * @param id Local Vertex which should be removed from the graph
    * @param f functor whose operator(GlobalEdge) is called for every removed edge
    **/
    template<typename Functor>
    void removeVertex(LocalVertex id, Functor& f);
    //no default template arguments for template functions allowed before c++0x, so a little workaround
    void removeVertex(LocalVertex id) ;

    /**
    * @brief Removes a vertex from the cluster or it's subclusters and applys functor to removed edges
    *
    * Removes the vertex from the graph or subclusters and invalidates the global vertex id. Also all edges connecting
    * to this vertex will be removed (upstream and downstream) after the functor was applied to them. The functor
    * needs to implement operato()(LocalEdge edge).
    *
    * @param id Global Vertex which should be removed from the graph
    * @param f functor whose operator(LocalEdge) is called on every removed edge
    **/
    template<typename Functor>
    void removeVertex(GlobalVertex id, Functor& f);
    //no default template arguments for template functions allowed before c++0x, so a little workaround
    void removeVertex(GlobalVertex id);

    /**
    * @brief Removes a global Edge from the cluster or it's subclusters
    *
    * Removes the edge from the graph or subclusters and invalidates the global edge id. If the local edge holds
    * only this global one it will be removed also.
    *
    * @param id Global Edge which should be removed from the graph
    * @return bool indicates if the global id could be removed
    **/
    void removeEdge(GlobalEdge id);

    /**
    * @brief Removes a local edge from the cluster and calls the functor for all removed global edges
    *
    * Removes the edge from the graph and invalidates the global edges. The Functor needs to provide
    * operator()(GlobalEdge). If no functor is needed just use boost::remove_edge.
    *
    * @param id Global Edge which should be removed from the graph
    * @param f functor whoms operator(GlobalEdge) is called
    * @return bool indicates if the global id could be removed
    **/
    template<typename Functor>
    void removeEdge(LocalEdge id, Functor& f);


    
    /********************************************************
    * Vertex and Cluster moving
    * *******************************************************/

    /**
     * @brief Move a vertex to a direct child subcluster
     *
     * This function moves the LocalVertex to the subcluster and reconnects all other vertices and clusters. The
     * moved vertex will hold it's global descriptor but get a new local one assigned (the one returned). The same
     * stands for all edges which use the moved vertex: global descriptors stay the same, but they are moved to new
     * local edges. It's allowed to move cluster vertices with this function.
     * @note: The subcluster must be a direct child of this cluster, it is not allowed to be further down in the 
     *        hirarchy. It must also not be a parent cluster.
     * @param v the LocalVertex to be moved
     * @param cg reference to the subcluster to which v should be moved
     * @return LocalVertex the local descriptor of the moved vertex in the subcluster it was moved to
     **/
    LocalVertex moveToSubcluster(LocalVertex v, std::shared_ptr<ClusterGraph> cg);

    /**
     * @brief Move a vertex to a direct child subcluster
     *
     * Overloaded convinience function which fetches the the cluster reference for the local descriptor and calls
     * the correct equivalent.
     * @note The subcluster defined by the local vertex must be a direct child of this cluster, no deeper 
     *       nested hirarchies are allowed
     * @param v the LocalVertex to be moved
     * @param Cluster the local vertex descriptor representing the subcluster to which v should be moved
     * @return LocalVertex the local descriptor of the moved vertex in the subcluster it was moved to
     **/
    LocalVertex moveToSubcluster(LocalVertex v, LocalVertex Cluster);

    /**
     * @brief Move a vertex to any child subcluster
     *
     * This function does move the vertex to a child exactly like the other two overloads, however, it 
     * is allowed to specify a global vertex which is deeper down in a hirarchy, meaning it must not 
     * be a direct child of this cluster. This allows to move local vertices over multiple subcluster 
     * borders with a single call. The global vertex is the vertex descriptor of the cluster in its own parent. 
     * For example, if you have cluster1 -> cluster2 ->cluster2 and want to move something from cluster1 
     * to cluster 3, the GlobalVertex is the cluster vertex within cluster 2.
     *
     * @param v the LocalVertex to be moved
     * @param Cluster the global vertex descriptor representing the subcluster to which v should be moved
     * @return LocalVertex the local descriptor of the moved vertex in the subcluster it was moved to
     **/
    LocalVertex moveToSubcluster(LocalVertex v, GlobalVertex Cluster);


    /**
     * @brief Move a vertex to the parent cluster.
     *
     * This function moves a vertex one step up in the subcluster hierarchie and reconnects all other vertices and clusters.
     * The moved vertex will hold it's global descriptor but get a new local one assigned (the one returned). The same
     * stands for all edges which use the moved vertex: global descriptors stay the same, but they are moved to new
     * local edges. Note that this function is the inverse of moveToSubcluster, and doing Pseudocode:
     * moveToParent(moveToSubcluster(v)) does nothing (only the local descriptor of the moved vertex is
     * diffrent afterwards).
     *
     * @param v Local vertex which should be moved to the parents cluster
     * @return LocalVertex the local descriptor of the moved vertex, valid in the parent cluster only.
     **/
    LocalVertex moveToParent(LocalVertex v);


    /********************************************************
    * Stuff
    * *******************************************************/

 private:
    ClusterMap  m_clusters;
    Graph m_graph;
    std::weak_ptr<ClusterGraph> m_parent;
    IDpointer 	  m_id;
    bool copy_mode; //no changing itself when copying

    /* Searches the global vertex in all local vertices of this graph, and returns the local
     * one which holds the global vertex. If not successfull the local vertex returned will be
     * invalid and the bool parameter will be false. If recursive = true, all subclusters will
     * be seached too, however, if found there the retourned local vertex will be the vertex
     * representing the toplevel cluster holding the global vertex in the initial graph.
     * */
    std::pair<LocalVertex, bool> getContainingVertex(GlobalVertex id, bool recursive = true);

    /* Searches the local vertex holding the specified global one in this and all it's subclusters.
     * If found, the holding local vertex and the graph in which it is valid will be returned.
     * */
    fusion::vector<LocalVertex, std::shared_ptr<ClusterGraph>, bool> getContainingVertexGraph(GlobalVertex id);
    
    /* Searches the global edge in all local edges of this graph, and returns the local
     * one which holds the global edge. If not successfull the local edge returned will be
     * invalid and the bool parameter will be false.
     * */
    std::pair<LocalEdge, bool> getContainingEdge(GlobalEdge id);
    
    fusion::vector<LocalEdge, std::shared_ptr<ClusterGraph>, bool> getContainingEdgeGraph(GlobalEdge id);

public:
    //may hold properties which have Eigen3 objects and therefore need alignment
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

/** @} */

//***************************************
//functors needed for implementation only
//***************************************

template<typename Functor, typename Graph>
struct apply_remove_prediacte {
    Functor& func;
    GlobalVertex vert;
    GlobalEdge edge;
    bool isEdge;

    apply_remove_prediacte(Functor& f, GlobalVertex v) : func(f), vert(v), isEdge(false) {};
    apply_remove_prediacte(Functor& f, GlobalEdge e) : func(f), edge(e), vert(0), isEdge(true) {};
    bool operator()(typename Graph::edge_bundle_single& e) {
        bool res;

        //this predicate can be used to compare the edge itself or the vertives it connects. See
        //if we are a relevant edge
        if(isEdge)
            res = (edge == e.template getProperty<EdgeProperty>());
        else
            res = (vert == e.template getProperty<EdgeProperty>().source) || (vert == e.template getProperty<EdgeProperty>().target);

        //we are a hit, invoke the functor.
        if(res || vert < 0)
            func(e.template getProperty<EdgeProperty>());

        return res || vert < 0;
    }
};

//***********************
//Function implementation
//***********************


template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
template<typename Functor>
void ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::copyInto(std::shared_ptr<ClusterGraph> into, Functor& functor) const {

    //lists does not provide vertex index, so we have to build our own (cant use the internal
    //vertex_index property as we would need to reset the indices and that's not possible in const graph)
    typedef std::map<LocalVertex, int> IndexMap;
    IndexMap mapIndex;
    boost::associative_property_map<IndexMap> propmapIndex(mapIndex);

    std::pair<local_vertex_iterator, local_vertex_iterator>  vit = boost::vertices(m_graph);

    for(int c = 0; vit.first != vit.second; vit.first++, c++)
        put(propmapIndex, *vit.first, c);

    //first copy all vertices and edges, but be aware that the objects in the new graph
    //are also copys only and point to the old graph. there is a bug in older boost version
    //(<1.5 i belive) that breaks vertex_all propety map for bundled properties, so we
    //have to create our own copie functors
    into->clear();
    vertex_copier<Graph> vc(m_graph, *into);
    edge_copier<Graph> ec(m_graph, *into);
    boost::copy_graph(m_graph, *into, boost::vertex_index_map(propmapIndex).vertex_copy(vc).edge_copy(ec));

    //set the IDgen to the same value to avoid duplicate id's in the copied cluster
    into->m_id->setCount(m_id->count());

    //now that we have all vertices we can recreate the subclusters
    std::pair<const_cluster_iterator, const_cluster_iterator> it = clusters();

    for(; it.first != it.second; it.first++) {
        //create the new Graph
        std::shared_ptr<ClusterGraph> ng = std::shared_ptr<ClusterGraph> (new ClusterGraph(into));

        //we already have the new vertex, however, we need to find it
        GlobalVertex gv = Base::getGlobalVertex((*it.first).first);
        LocalVertex  lv = into->getLocalVertex(gv).first;

        //add the new graph to the subclustermap
        into->m_clusters[lv] = ng;

        //copy the subcluster
        (*it.first).second->copyInto(ng, functor);
    }

    //lets see if the objects need special treatment
    into->for_each_object(functor, false);
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
void ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::setCopyMode(bool on) {
    copy_mode = on;
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
template<typename P>
typename P::type& ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::getSubclusterProperty(LocalVertex v) {
    return getVertexCluster(v)->template getProperty<P>();
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
std::pair<std::shared_ptr< ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop> >, LocalVertex> ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::createCluster() {
    typename Base::vertex_bundle vp;
    vp.template setProperty<VertexProperty>(m_id->generate());
    LocalVertex v = boost::add_vertex(vp, m_graph);
    std::shared_ptr<ClusterGraph> sp = std::static_pointer_cast<ClusterGraph>(sp_base::shared_from_this());
    return std::pair<std::shared_ptr<ClusterGraph>, LocalVertex> (m_clusters[v] = std::shared_ptr<ClusterGraph> (new ClusterGraph(sp)), v);
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
inline std::shared_ptr< ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop> >
ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>:: parent()  {
    return std::shared_ptr<ClusterGraph> (m_parent);
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
inline const std::shared_ptr< ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop> >
ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::parent() const     {
    return std::shared_ptr<ClusterGraph> (m_parent);
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
bool ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::isRoot() const {
    return m_parent.expired();
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
std::shared_ptr< ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop> >
ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::root()             {

    return isRoot() ? std::static_pointer_cast<ClusterGraph>(sp_base::shared_from_this()) : parent()->root();
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
const std::shared_ptr< ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop> >
ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::root() const    {
    return isRoot() ? sp_base::shared_from_this() : parent()->root();
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
std::pair<typename ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::cluster_iterator, typename ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::cluster_iterator>
ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::clusters() {
    return std::make_pair(m_clusters.begin(), m_clusters.end());
}

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
std::pair<typename ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::const_cluster_iterator, typename ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::const_cluster_iterator>
ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::clusters() const {
    return std::make_pair(m_clusters.begin(), m_clusters.end());
}

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
std::size_t
ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::numClusters() const {
    return m_clusters.size();
}

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
bool ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::isCluster(const LocalVertex v) const {
    return (m_clusters.find(v) != m_clusters.end());
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
std::shared_ptr< ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop> >
ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::getVertexCluster(LocalVertex v) {
    if(isCluster(v))
        return m_clusters[v];

    //TODO:throw if not a cluster
    return std::shared_ptr< ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop> >();//sp_base::shared_from_this();
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
LocalVertex     ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::getClusterVertex(std::shared_ptr<ClusterGraph> g) {
    std::pair<cluster_iterator, cluster_iterator> it = clusters();

    for(; it.first != it.second; it.first++) {
        if((*it.first).second == g)
            return (*it.first).first;
    }

    throw graph::cluster_error() <<  boost::errinfo_errno(12) << error_message("Cluster is not part of this graph");
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
template<typename Functor>
void ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::removeCluster(std::shared_ptr<ClusterGraph> g, Functor& f) {
    removeCluster(getClusterVertex(g), f);
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
void ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::removeCluster(std::shared_ptr<ClusterGraph> g) {
    placehoder p;
    removeCluster(getClusterVertex(g), p);
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
void ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::clearClusters() {
    m_clusters.clear();
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
template<typename Functor>
void ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::removeCluster(LocalVertex v, Functor& f) {

    typename ClusterMap::iterator it = m_clusters.find(v);

    if(it == m_clusters.end())
        throw graph::cluster_error() <<  boost::errinfo_errno(11) << error_message("Cluster is not part of this graph");

    std::pair<LocalVertex, std::shared_ptr<ClusterGraph> > res = *it;

    //apply functor to all vertices and edges in the subclusters
    f(res.second);
    res.second->remove_vertices(f, true);

    //remove from map, delete subcluster and remove vertex
    m_clusters.erase(v);
    boost::clear_vertex(v, m_graph);    //should not be needed, just to ensure it
    boost::remove_vertex(v, m_graph);
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
void ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::removeCluster(LocalVertex v) {
    placehoder p;
    removeCluster(v, p);
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
template<typename Functor>
void ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::remove_vertices(Functor& f, bool recursive) {

    std::pair<local_vertex_iterator, local_vertex_iterator>  vit = boost::vertices(m_graph);

    //we iterate forward before deleting to not invalidate our iterator
    while(vit.first != vit.second) {
        LocalVertex v = * (vit.first);
        vit.first++;

        if(!isCluster(v)) {
            //let the functor know we remove this vertex
            f(Base::getGlobalVertex(v));
            //need to do this to allow the removal of all relevant edges to this vertex, even upstream
            removeVertex(v, f);
        }
    };

    if(recursive) {
        cluster_iterator cit;

        for(cit = m_clusters.begin(); cit != m_clusters.end(); cit++) {
            f((*cit).second);
            (*cit).second->remove_vertices(f, recursive);
        }
    }
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
fusion::vector<LocalVertex, GlobalVertex>
ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::addVertex() {

    vertex_bundle vp;
    vp.template setProperty<VertexProperty>(m_id->generate());
    LocalVertex v = boost::add_vertex(vp, m_graph);

    return fusion::make_vector(v, m_id->count());
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
fusion::vector<LocalVertex, GlobalVertex>
ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::addVertex(GlobalVertex gv) {

    std::pair<LocalVertex, bool> res = Base::getLocalVertex(gv);

    if(!res.second) {
        vertex_bundle vp;
        vp.template setProperty<VertexProperty>(gv);
        LocalVertex v = boost::add_vertex(vp, m_graph);

        //ensure that we never create this id, as it is used now
        if(gv > m_id->count())
            m_id->setCount(gv);

        return fusion::make_vector(v, gv);
    };

    return fusion::make_vector(res.first, gv);
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
fusion::vector<LocalEdge, GlobalEdge, bool>
ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::addEdge(LocalVertex source, LocalVertex target) {

    //manual edge creation with cluster is not allowed
    if((source == target) || isCluster(source) || isCluster(target))
        return fusion::make_vector(LocalEdge(), GlobalEdge(), false);

    LocalEdge e;
    bool done;
    boost::tie(e, done) = boost::edge(source, target, m_graph);

    //if done=true the edge alredy existed
    if(!done)
        boost::tie(e, done) = boost::add_edge(source, target, m_graph);

    if(!done)
        return fusion::make_vector(LocalEdge(), GlobalEdge(), false);

    //init the bundle corecctly for new edge
    GlobalEdge global = { m_graph[source].template getProperty<VertexProperty>(),
                          m_graph[target].template getProperty<VertexProperty>(), m_id->generate() };
    edge_bundle_single s;
    s.template setProperty<EdgeProperty>(global);
    auto& vec = m_graph[e].template getPropertyAccessible<GEdgeProperty>();
    vec.push_back(s);
    m_graph[e].template markPropertyChanged<GEdgeProperty>();

    return fusion::make_vector(e, global, true);
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
fusion::vector<LocalEdge, GlobalEdge, bool, bool>
ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::addEdge(GlobalVertex source, GlobalVertex target) {

    LocalVertex v1, v2;
    LocalEdge e;
    bool d1, d2, d3;
    boost::tie(v1, d1) = getContainingVertex(source);
    boost::tie(v2, d2) = getContainingVertex(target);

    //if one vertex is not accessible from here this function fails
    if(!(d1 && d2))
        return fusion::make_vector(LocalEdge(), GlobalEdge(), false, false);

    //if both vertices are in a subcluster this one must do the job as we cant access the local edge from here
    if(v1 == v2 && isCluster(v1)) {
        fusion::vector<LocalEdge, GlobalEdge, bool, bool> res = getVertexCluster(v1)->addEdge(source, target);
        fusion::at_c<3> (res) = false;
        return res;
    }

    //check if we already have that Local edge
    boost::tie(e, d3) = boost::edge(v1, v2, m_graph);

    if(!d3)
        boost::tie(e, d3) = boost::add_edge(v1, v2, m_graph);

    if(!d3)
        return fusion::make_vector(LocalEdge(), GlobalEdge(), false, false);

    //init the bundle corectly for new edge
    GlobalEdge global = { source, target, m_id->generate() };
    edge_bundle_single s;
    s.template setProperty<EdgeProperty>(global);
    auto& vec = m_graph[e].template getPropertyAccessible<GEdgeProperty>();
    vec.push_back(s);
    m_graph[e].template markPropertyChanged<GEdgeProperty>();
    
    return fusion::make_vector(e, global, true, true);

};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
fusion::vector<LocalEdge, GlobalEdge, bool, bool>
ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::addEdgeGlobal(GlobalVertex source, GlobalVertex target) {
    return addEdge(source, target);
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
fusion::vector<LocalEdge, ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>*, bool>
ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::getLocalEdgeGraph(GlobalEdge e) {
    return getContainingEdgeGraph(e);
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
fusion::vector<LocalVertex, std::shared_ptr< ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop> >, bool>
ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::getLocalVertexGraph(GlobalVertex v) {
    return getContainingVertexGraph(v);
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
template<typename Functor>
void ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::downstreamRemoveVertex(GlobalVertex v, Functor& f) {

    std::pair<LocalVertex, bool> res = getContainingVertex(v);

    //we don't throw, as this function gets invoked recursivly and it may happen that the
    //vertex to remove is only in the top layers, not the button ones
    if(!res.second)
        return;


    //iterate over every edge that connects to the global vertex or the cluster in which it is in
    std::vector<LocalEdge> re; //remove edges
    std::pair<local_out_edge_iterator,  local_out_edge_iterator> it = boost::out_edges(res.first, m_graph);

    for(; it.first != it.second; it.first++) {
        auto& vec = m_graph[* (it.first)].template getPropertyAccessible<GEdgeProperty>();
        vec.erase(std::remove_if(vec.begin(), vec.end(), apply_remove_prediacte<Functor, ClusterGraph> (f, v)), vec.end());
        m_graph[* (it.first)].template markPropertyChanged<GEdgeProperty>();
        
        if(vec.empty())
            re.push_back(* (it.first));
    };

    std::for_each(re.begin(), re.end(), boost::bind(&ClusterGraph::simpleRemoveEdge, this, _1));

    //if we have the real vertex here and not only a containing cluster we can delete it
    if(!isCluster(res.first)) {
        boost::clear_vertex(res.first, m_graph);    //just to make sure, should be done already
        boost::remove_vertex(res.first, m_graph);
    };

    //lets go downstream
    for(cluster_iterator it = m_clusters.begin(); it != m_clusters.end(); it++)
        ((*it).second)->downstreamRemoveVertex(v, f);
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
void ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::simpleRemoveEdge(LocalEdge e) {
    boost::remove_edge(e, m_graph);
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
template<typename Functor>
void ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::removeVertex(LocalVertex id, Functor& f) {
    //it is important to delete the global vertex, not the only local one as it's possible that
    //we are in a subcluster and there are connections to the global vertex in the parent. They
    //need to be deleted too.
    removeVertex(Base::getGlobalVertex(id), f);
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
void ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::removeVertex(LocalVertex id) {
    placehoder p;
    removeVertex(Base::getGlobalVertex(id), p);
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
template<typename Functor>
void ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::removeVertex(GlobalVertex id, Functor& f) {
    root()->downstreamRemoveVertex(id, f);
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
void ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::removeVertex(GlobalVertex id) {
    placehoder p;
    removeVertex(id, p);
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
void ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::removeEdge(GlobalEdge id) {
    
    fusion::vector<LocalEdge, std::shared_ptr<ClusterGraph>, bool> res = getContainingEdgeGraph(id);

    if(!fusion::at_c<2> (res))
        return; //TODO:throw

    placehoder p;
    auto& vec = ((*fusion::at_c<1> (res)) [fusion::at_c<0> (res)]).template getPropertyAccessible<GEdgeProperty>();
    vec.erase(std::remove_if(vec.begin(), vec.end(), apply_remove_prediacte<placehoder, ClusterGraph> (p, id)), vec.end());
    ((*fusion::at_c<1> (res)) [fusion::at_c<0> (res)]).template markPropertyChanged<GEdgeProperty>();
    
    if(vec.empty())
        boost::remove_edge(fusion::at_c<0> (res), fusion::at_c<1> (res)->getDirectAccess());    
    
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
template<typename Functor>
void ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::removeEdge(LocalEdge id, Functor& f) {

    auto& vec = m_graph[id].template getPropertyAccessible<GEdgeProperty>();
    std::for_each(vec.begin(), vec.end(), boost::bind<void> (boost::ref(apply_remove_prediacte<placehoder, ClusterGraph> (f, -1)), _1));
    m_graph[id].template markPropertyChanged<GEdgeProperty>();
    boost::remove_edge(id, m_graph);
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
LocalVertex
ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::moveToSubcluster(LocalVertex v, GlobalVertex gv) {

    auto res = getContainingVertex(gv); //std::pair<LocalVertex, bool> 
    dcm_assert(res.second);
    //we definitely move it to the subcluster
    auto nv = moveToSubcluster(v, res.first);
    
    //now check if this is already the correct one, and if not move it further down.
    if(Base::getGlobalVertex(res.first) != gv)        
        return getVertexCluster(res.first)->moveToSubcluster(nv, gv);
    
    return nv;
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
LocalVertex
ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::moveToSubcluster(LocalVertex v, LocalVertex Cluster) {

    std::shared_ptr<ClusterGraph> cg = getVertexCluster(Cluster);
    return moveToSubcluster(v, cg);
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
LocalVertex
ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::moveToSubcluster(LocalVertex v, std::shared_ptr<ClusterGraph> cg) {

    LocalVertex Cluster = getClusterVertex(cg);
    
    std::pair<local_out_edge_iterator, local_out_edge_iterator> it =  boost::out_edges(v, m_graph);

    /* add the later removed edges to the coressponding existing edges
     * (or create new edges between adjacent vertices of moved vertex and cluster).
     * also get the edge between cluster and vertex while iterating */
    for(; it.first != it.second; it.first++) {

        LocalVertex target = boost::target(*it.first, m_graph);
        if(target == v )
            target = boost::source(*it.first, m_graph);

        if(target != Cluster) {

            //get or create the edge between the old edge target and the cluster
            LocalEdge e;
            bool done;
            boost::tie(e, done) = boost::edge(target, Cluster, m_graph);

            if(!done)
                boost::tie(e, done) = boost::add_edge(target, Cluster, m_graph);

            //if(!done) TODO: throw

            auto& ep = m_graph[*it.first].template getPropertyAccessible<GEdgeProperty>();
            auto& nep = m_graph[e].template getPropertyAccessible<GEdgeProperty>();
            nep.insert(nep.end(), ep.begin(), ep.end());
            m_graph[e].template markPropertyChanged<GEdgeProperty>();
        }
    }

    /* Create new Vertex in Cluster and map the edge to vertices and clusters in the cluster
    * if a connection existed */
    LocalVertex nv = boost::add_vertex(m_graph[v], cg->getDirectAccess());

    //resort cluster parentship if needed
    if(isCluster(v)) {

        cg->m_clusters[nv] = m_clusters[v];
        cg->m_clusters[nv]->m_parent = cg;
        m_clusters.erase(v);
    }

    std::pair<LocalEdge, bool> moveedge = boost::edge(v, Cluster, m_graph);

    if(moveedge.second) {
        auto& vec = m_graph[moveedge.first].template getPropertyAccessible<GEdgeProperty>();

        for(edge_single_iterator i = vec.begin(); i != vec.end(); i++) {

            //get the global vertex to which the global edge points and find the local vertex holding this
            //global one
            GlobalEdge global = typename Base::global_extractor()(*i);
            GlobalVertex target;
            //bit cumbersome to support moving clusters
            target = (cg->getContainingVertex(global.source).first == nv) ? global.target : global.source;
            std::pair<LocalVertex, bool> res = cg->getContainingVertex(target);
            //if(!res.second) TODO: throw

            //get or create the edge between the new vertex and the target
            LocalEdge e;
            bool done;
            boost::tie(e, done) = boost::edge(nv, res.first, cg->getDirectAccess());

            if(!done)
                boost::tie(e, done) = boost::add_edge(nv, res.first, cg->getDirectAccess());

            //if(!done) TODO: throw

            //push the global edge to the local edge
            auto& gvec = (*cg)[e].template getPropertyAccessible<GEdgeProperty>();
            gvec.push_back(*i);
            (*cg)[e].template markPropertyChanged<GEdgeProperty>();
        };
    }

    //all global edges concerning the move vertex are processed and it is moved to the subcluster,
    //lets destroy it in the local cluster
    boost::clear_vertex(v, m_graph);
    boost::remove_vertex(v, m_graph);

    return nv;
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
LocalVertex
ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::moveToParent(LocalVertex v) {

    //if(isRoot()) TODO:throw

    //create new vertex
    vertex_bundle& vb = m_graph[v];
    LocalVertex nv = boost::add_vertex(vb, parent()->getDirectAccess());

    //regrouping if needed
    if(isCluster(v)) {
        parent()->m_clusters[nv] = m_clusters[v];
        parent()->m_clusters[nv]->m_parent = m_parent;
        m_clusters.erase(v);
    }

    GlobalVertex gv = vb.template getProperty<VertexProperty>();

    //get all out_edges of this cluster in the parentcluster (because only they can hold relevant global_Edgs)
    std::vector<LocalEdge> edge_vec;
    std::shared_ptr<ClusterGraph> sp = std::static_pointer_cast<ClusterGraph>(sp_base::shared_from_this());
    LocalVertex this_v = parent()->getClusterVertex(sp);
    std::pair<local_out_edge_iterator, local_out_edge_iterator> it = boost::out_edges(this_v, parent()->getDirectAccess());

    for(; it.first != it.second; it.first++) {
        //iterate all global edges and find relevant ones
        auto& vec = ((parent()->getDirectAccess()) [*it.first]).template getPropertyAccessible<GEdgeProperty>();
        edge_single_iterator i = vec.begin();

        while(i != vec.end()) {

            GlobalEdge global = typename Base::global_extractor()(*i);
            GlobalVertex target;

            //a bit cumbersome to allow cluster moving
            if(parent()->getContainingVertex(global.source).first == nv)
                target = global.target;
            else if(parent()->getContainingVertex(global.target).first == nv)
                target = global.source;
            else {
                i++;
                continue;
            }

            std::pair<LocalVertex, bool> res = parent()->getContainingVertex(target);

            //get or create the edge between the new vertex and the target
            LocalEdge e;
            bool done;
            boost::tie(e, done) = boost::edge(nv, res.first, parent()->getDirectAccess());

            if(!done)
                boost::tie(e, done) = boost::add_edge(nv, res.first, parent()->getDirectAccess());

            //if(!done) TODO: throw

            //push the global edge bundle to the new local edge and erase it in the old
            auto& gvec =  (parent()->getDirectAccess())[e].template getPropertyAccessible<GEdgeProperty>();
            gvec.push_back(*i);
            (parent()->getDirectAccess())[e].template markPropertyChanged<GEdgeProperty>();
            
            i = vec.erase(i);
        }

        //see if we should destroy this edge (no global edges remain in local one)
        if(vec.empty())
            edge_vec.push_back(*it.first);
    }

    //create a edge between new vertex and this cluster and add all global edges from within this cluster
    it = boost::out_edges(v, m_graph);
    LocalEdge e;

    if(it.first != it.second)
        e = boost::add_edge(nv, this_v, parent()->getDirectAccess()).first;

    for(; it.first != it.second; it.first++) {
        auto& ep  = m_graph[*it.first].template getPropertyAccessible<GEdgeProperty>();
        auto& nep = (*parent())[e].template getPropertyAccessible<GEdgeProperty>();
        nep.insert(nep.end(), ep.begin(), ep.end());
        
        (*parent())[e].template markPropertyChanged<GEdgeProperty>();
    }

    //all global edges concerning the move vertex are processed and it is moved to the parent,
    //lets destroy it in the local cluster
    boost::clear_vertex(v, m_graph);
    boost::remove_vertex(v, m_graph);

    //it's possible that some local edges in the parent are empty now, let's destroy them
    for(std::vector<LocalEdge>::iterator it = edge_vec.begin(); it != edge_vec.end(); it++)
        boost::remove_edge(*it, parent()->getDirectAccess());

    return nv;
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
std::pair<LocalVertex, bool>
ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::getContainingVertex(GlobalVertex id, bool recursive) {

    //check all vertices if they are the id
    std::pair<local_vertex_iterator, local_vertex_iterator>  it = boost::vertices(m_graph);

    for(; it.first != it.second; it.first++) {
        if(id == m_graph[*it.first].template getProperty<VertexProperty>())
            return std::make_pair(*it.first, true);
    }

    //check all clusters if they have the id
    if(recursive) {
        for(cluster_iterator it = m_clusters.begin(); it != m_clusters.end(); it++) {
            std::pair<LocalVertex, bool> res = ((*it).second)->getContainingVertex(id);

            if(res.second)
                return std::make_pair((*it).first, true);
        }
    }

    return std::make_pair((LocalVertex) NULL, false);
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
fusion::vector<LocalVertex, std::shared_ptr< ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop> >, bool>
ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::getContainingVertexGraph(GlobalVertex id) {

    LocalVertex v;
    bool done;
    boost::tie(v, done) = getContainingVertex(id);

    if(!done)
        return fusion::make_vector(LocalVertex(), std::shared_ptr<ClusterGraph>(), false);

    if(isCluster(v) && (Base::getGlobalVertex(v) != id))
        return m_clusters[v]->getContainingVertexGraph(id);
    else
        return fusion::make_vector(v, std::static_pointer_cast<ClusterGraph>(sp_base::shared_from_this()), true);
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
std::pair<LocalEdge, bool>
ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::getContainingEdge(GlobalEdge id) {

    LocalVertex v1, v2;
    bool d1, d2;
    boost::tie(v1, d1) = getContainingVertex(id.source, true);
    boost::tie(v2, d2) = getContainingVertex(id.target, true);

    if(!((d1 && d2) && (v1 != v2)))
        return std::make_pair(LocalEdge(), false);

    return boost::edge(v1, v2, m_graph);
};

template< typename edge_prop, typename globaledge_prop, typename vertex_prop, typename cluster_prop>
fusion::vector<LocalEdge, std::shared_ptr< ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop> >, bool>
ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop>::getContainingEdgeGraph(GlobalEdge id) {

    LocalVertex v1, v2;
    bool d1, d2;
    boost::tie(v1, d1) = getContainingVertex(id.source, true);
    boost::tie(v2, d2) = getContainingVertex(id.target, true);

    if(!(d1 && d2))
        return fusion::make_vector(LocalEdge(), std::shared_ptr< ClusterGraph<edge_prop, globaledge_prop, vertex_prop, cluster_prop> >(), false);

    if(v1 == v2)
        return m_clusters[v1]->getContainingEdgeGraph(id);

    std::shared_ptr<ClusterGraph> sp = std::static_pointer_cast<ClusterGraph>(sp_base::shared_from_this());
    return fusion::make_vector(boost::edge(v1, v2, m_graph).first, sp, true);
};

} //namespace graph
} //namespace dcm


#endif // CLUSTERGRAPH_HPP





