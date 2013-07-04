/*
    openDCM, dimensional constraint manager
    Copyright (C) 2013  Stefan Troeger <stefantroeger@gmx.net>

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

/** @addtogroup Core
 * @{*/

/** @defgroup ClusterGraph Custom Graph
 *
 * @brief A custom graph class based on boost::adjacency_list with support for 'clustering' and property
 * and object support.
 * 
 * \section graph Intruducing the ClusterGraph
 * 
 * To allow advanced processing of a given geometric constellation a modell is needed which holds all the
 * information in a accessible manner and allows fast and complex interaction. A graph is perfectöy suited 
 * for such a task and therefore used to represent the users geometric problem.
 * 
 * In the scope of this library vertices are used as geometry and edges as the constraint's between them.
 * This gives an abstract representation of the given dimensional constraint problem and alows further 
 * analysis. However, in this context we quickly reach a fundemental limitation: We cannot group a small
 * part of the graph in a non-disruptive way. The need for grouping arises in the fact, that fully constraint
 * geometrys can be handled as rigid system with a rotation and translation. Thats very helpfull for reducing
 * the constraint system in it's complexity. With a boost graph however, one would need to mark all edges and 
 * vertices as beeing part of a group and handle them in all algorithms seperatly. Thats a huge inconvienience.
 * An alternative is the boost subgraph code. This allows to group parts of a graph and use it as subgraph, 
 * apearing in the original graph as a single vertex. The problem that arises here is that all edges between
 * the vertices of the parent graph to the vertices in the subgraph get lost and therefore the given constraint
 * system would not be represented anymore. This problem is solved by this class too.
 *
 * The ClusterGraph extends the boost graph in two ways:
 * First it allows to specify properties and objects which are stored at vertices, edges and the
 * graph itsef. There are functions to compfortable access and set them.
 * Second, and most important, it allows subclustering. That means creating independend graphs as
 * child, where every vertex can be accessed by edges from the parent graph. A ClusterGraph is noncopyable, it
 * is therefore always passed around as shared_ptr.
 * 
 * \subsection clustergraph The cluster concept
 * 
 * The cluster graph answers to the grouping problem by allowing to 'cluster' it. A cluster is a part of the
 * whole graph, that can be treated independently. Imagine this graph: 
 * \image html graph.svg A simple graph
 * It has three vertices and 4 edges, vertex 1 and 2 are connected by two edges. Now imagine we want to group vertex 2 and
 * 3 and represent them as single vertex in our graph without loosing the edges. That would inevitable result
 * in creating a new graph where we put vertex 2 and 3 as well as edge 3 in. But here it goes: a new graph
 * means new vertex and edge descriptors and it would not be possible to access our entitys in the new graph 
 * with the old descriptors. Therefoe we first introduce a new descriptor system, a global one. Every edge and
 * every vertex has besides it's local boost descriptor a global one which indentifies it over all it's lifetime.
 * That would result in this new graph: 
 * \image html clustergraph.svg ClusterGraph with global descriptors
 * You can see how a global descriptor was added to every entity, but so see something else: Some strange global
 * edges appeared and are somewhat grouped in local edges. Why is that? Now thats because we need some way of 
 * collecting edges in a simple manner. In the original graph, we had two edges from vertex 1 to 2. In the cluster
 * graph there is only one real edge, but it holds the two connections as two global edges. The global edge is a
 * own type and holds the global edge descriptor, so don't confuse it with a descriptor itself. The hole purpose of
 * this construct will be evident later. 
 * 
 * At this point we have estalblished a way of addresssing our vertices and edges in a global and never-ever-changing 
 * manner. So no matter what for example the boost algorithms do to our entitys, the descriptors stay unchanged! (Note 
 * that the boost local descriptors may will stay unchanged as long as you do not move the entity out of the graph,
 * thats a result of the listS usage and therefore the local descriptors beeing viod pointers. That would be diffrent
 * for vector storage.) So with this in place, we can start our cluster process. That means we take vertex 1, vertex 2
 * and edge 3 and put them into a new ClusterGraph. Of course we still want to hold the new clustergraph in our current 
 * graph, as all of its entitys are part of our system. So we include it as a vertex of its own. 
 * 
 * When doing that we need to be careful: Whats about edges 1 and 2? After the change they would connect vertices in 
 * diffrent graphs? And if the new graph is a vertex in the original one now, how could the edges connect to the original
 * vertices? Simple: we have a local edge between vertex 1 and the new cluster vertex, so we can stack all global edges 
 * inside the local one! Not that obvious? Take a look here:
 * \image html subclustergraph.svg Clustergraph with vertex 1 and 2 in a subcluster
 * Thats the single reason we introduced the stacked local edges in the first place, to allow a arbitrary number of
 * global ones be stored in such a way, that a local analysy only sees one local edge. So in the local point of view,
 * the shown graph has two vertices and one edge connecting it.Thats what the boost algorithms see. In the global 
 * perspective however, all our 3 vertices and 4 edges are still there. And you could also stack it deeper, for example
 * creating a new subcluster in the already created one and therefore stacking 3 ClusterGraph's. You will still be able
 * to connect vertex 1 to all vertices in the new nested subclusterwithglobal edges, there is no restriction. And thats
 * how the cluster concept works.
 * 
 * \subsection properties Adding propeties and objects
 * 
 * It is possible to at data to the ClusterGraph's entitys by attatching multiple \ref Property 's to them. This includes 
 * local edges, local vertices and clusters. The vertex and edge properties are intended to be used in boost graph algorithms and 
 * shall be used in combination with the \ref property_map interface. They allow to store algorithmic information  at a
 * convienient place for later evaluation. Global edges don't have properties as they are not used in boost 
 * algorithms. It's also possible to add properties to a hole clustergraph. This is needed to give the subclusters
 * a meaning of some kind. With attached properties it's possible to differentiate diffrent use cases and distuingish
 * the relevant cluster.
 * 
 * To add properties to the three possible places they need to be passed as mpl::vectors in the appropriate template
 * parameter of the ClusterGraph. The Property kind is ignored by the graph and assignment to edges, vertices or cluster
 * is only based on the template parameter order. However, it's important that the parameters are mpl::vector's and that
 * every type in it is a property as describet in \ref Property .
 * 
 * 
 * 
 *@}/

/**@}*/

