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

#ifndef DCM_MODULE_3D_H
#define DCM_MODULE_3D_H

#include "opendcm/core/object.hpp"
#include "opendcm/core/module.hpp"
#include "opendcm/core/utilities.hpp"
#include "opendcm/core/geometry.hpp"
#include "opendcm/core/constraint.hpp"
#include "opendcm/core/typeadaption.hpp"

#include "geometry.hpp"
#include "reduction.hpp"
#include "cluster.hpp"

#include <type_traits>

struct symbol;
namespace dcm {
       
template<typename ... types>
struct Module3D {

    typedef boost::mpl::int_<5> ID;

    template<typename Final, typename Stacked>
    struct type : public Stacked {

        typedef typename Stacked::Kernel Kernel;
      
        type() : Stacked() {};        
        ~type() {}
        
        //we setup the reduction graph after the constructors are done, as only then the graph 
        //is available (as well as the virtual function to get it)
        void init() {
            
            //that is our responsibility in the init stack
            Stacked::init();
            
            //the default cluster reduction is a numeric::Geometry. This is not good, we want it to 
            //be our parametric geometry
            auto reduction = this->getReductionGraph(geometry::Cluster3<Kernel>::id(), geometry::Cluster3<Kernel>::id());
            reduction->replaceSourceNode(std::make_shared<reduction::NumericGeometryNode<Kernel, numeric::Cluster3<Kernel>>>());
            //TODO:  Do this for all graphs that end in a cluster! 
            
            //we create our default edge reduction trees
            int pointID = geometry::Point3<Kernel>::id();
            module3d::setupPointPointReduction<Kernel>(this->getReductionGraph(pointID, pointID));
        };
        
        /**
         * @brief Container for 3D user geometry
         * 
         * You can access the once set geometry in two ways. For one simply by using the get<>() function. However, this requires 
         * to know the type of the geometry that is currently hold which must be passed as template argument. You can inquery if 
         * a certain type is hold by using the holdsGeometryType functions. This could look like
         * @code
         * if(geometry->holdsGeometryType<MaPointClass>()) {
         *     auto myPoint = geometry->get<MyPointClass>();
         *     auto alternative = dcm::get<MyPointClass>(geometry);
         * @endcode
         * Alternatively it is possible to apply a visitor object which gets called with the stored geometry type. The visitor can 
         * have a operator() for each expected type, or even a template operator if all should be handled equally. Note that if 
         * you pass in the visitor as temporary object the operator() as well as the passed type must be defined const. 
         * Note also that the return type of the visitor must be defined via the dcm::visitor<> template.
         * @code 
         * struct Visitor1 : dcm::visitor<int> {
         *      int operator()(MyPointClass& geometry) {return 1;};
         * }
         * struct Visitor2 : dcm::visitor<int> {
         *      template<typename T>
         *      int operator()(const T& geometry) const {return 1;};
         * }
         * 
         * Visitor1 visitor;
         * int num = geometry->apply(visitor);
         * int num = geometry->apply(Visitor2()); //works as the operator() is const
         * @endcode
         */
        struct Geometry3D : public Stacked::Object, public utilities::Variant<types...> {

            typedef utilities::Variant<types...> InheritedV;
            typedef typename Stacked::Object     InheritedO;
                       
            typedef symbolic::GeometryProperty GeometryProperty;
            typedef graph::VertexProperty      VertexProperty;
            
            DCM_OBJECT_ADD_PROPERTIES( InheritedO, (VertexProperty) )
            
        public:            
            Geometry3D(Final* system) 
                : Stacked::Object( Final::template objectTypeID<typename Final::Geometry3D>::ID::value ),
                m_system(system), m_type(-1) {};
            
            /**
             * @brief Set the content and stores the user type
             * 
             * This function is used to set the content of the Geometry3D container. It extract the contents of 
             * the given \a geometry and stores it in an accessible manner for the solver. Furhtermore it stores
             * the supplied user type for later access. The exact behaviour of this feature depends on the specified
             * qualifiers of the storable types as given as the Module3D template arguments. If the type was given
             * without any special qualifiers then a copy of the given user type is strored.
            */
            template<typename T>
            void set(const T& geometry) {
                
                typedef typename Final::Kernel  Kernel;
                typedef typename Kernel::Scalar Scalar;
                typedef geometry::extractor<typename geometry_traits<T>::type> extractor;
                typedef typename extractor::template primitive<typename Final::Kernel> Geometry;
                typedef symbolic::TypeGeometry<Geometry> TypeGeometry;

                BOOST_MPL_ASSERT((mpl::contains<mpl::vector<types...>, T>));
                
                //we may need to setup the graph. We do this here to allow to clear the geometry and later reinitialize by setting a new geometry.
                std::shared_ptr<typename Final::Graph> cluster = std::static_pointer_cast<typename Final::Graph>(m_system->getGraph());
                if(!holdsGeometry()) {                    
                    fusion::vector<graph::LocalVertex, graph::GlobalVertex> res = cluster->addVertex();
                    cluster->template setProperty<details::GraphObjectProperty>(fusion::at_c<0>(res), InheritedO::shared_from_this());
                    setVertexProperty(fusion::at_c<1>(res));
                };
                
                //store the type
                m_type = Geometry::id();
                InheritedV::m_variant = geometry;
                
            };
                        
            /**
             * @brief Check if the Geometry3D holds a geometry type
             * 
             * As a Geometry3D is only a container for user geometry it is only valid if a user type was set. 
             * To check if this already happend this function can be used. If it returns \a false then the 
             * \ref set function neets to be used to initialize the container before using it to setup a system.
             * \note When the geometry is created through the systems \ref addGeometry3D funcion it is always 
             * propperly initialized, false is only possible if the user creates the Geometry3D itself.
             * 
             * @return bool true if a valid geometry type was set
             */
            bool holdsGeometry() {
                return InheritedV::holdsType();
            };
            
            /**
             * @brief Check if a certain geometry type is stored
             * 
             * As a Geometry3D is only a container for user types it is often needed to query the type currently 
             * stored. This can be done with this function in two ways. First it is possible to query by the dcm
             * geometry type, for example dcm::Point3. Valid types are the one used in the geometry traits. The 
             * alternative is to directly provide the user types. Valid values are the ones used to construct the
             * Module3D. See the following example for the usage (given MyCustomPointType was regirstered before 
             * as dcm::Point3 via geometry_traits)
             * @code
             * boost::shared_ptr<Geometry3D> g(new Geometry3D(system));
             * g->set(MyCustomPointType())
             * assert(g->holdsGeometryType<dcm::Point3>())
             * assert(!g->holdsGeometryType<dcm::Line3>())
             * assert(g->holdsGeometryType<MyCustomPointType>())
             * assert(!g->holdsGeometryType<MyCustomLineType>())
             * @endcode
             * \note Unregisterd user types are not supported, the function cannot be compiled in such a 
             * situation
             */
            template<typename T>
            typename boost::enable_if<mpl::contains<mpl::vector<types...>, T>, bool>::type 
            holdsGeometryType() {
                return m_type == geometry::extractor<typename geometry_traits<T>::type>::template primitive<Kernel>::id();
            };
            template<typename T>
            typename boost::disable_if<mpl::contains<mpl::vector<types...>, T>, bool>::type 
            holdsGeometryType() {
                return m_type == geometry::extractor<T>::template primitive<Kernel>::id();
            };
           
            
        protected:
            Final* m_system;
            int    m_type;

            virtual void preprocessVertex(std::shared_ptr<graph::AccessGraphBase> g, 
                                          graph::LocalVertex lv, graph::GlobalVertex gv) override {
                
                auto cluster = std::static_pointer_cast<typename Final::Graph>(g);
                auto prop = cluster->template getProperty<GeometryProperty>(lv);
                
                //if the property is not available or does not hold the correct type we need to delete it, no way around
                if(prop && prop->getType() != m_type) {
                    delete prop;
                    prop = nullptr;
                }
                
                //if no property in existance we need to create it
                if(!prop) {
                    prop = InheritedV::apply(PropCreator());
                    cluster->template setProperty<GeometryProperty>(lv, prop);
                }
                    
                //we definitly need to set the value
                InheritedV::apply(PropAssigner(prop));
            };
            
            virtual void postprocessVertex(std::shared_ptr<graph::AccessGraphBase> g,
                                           graph::LocalVertex lv, graph::GlobalVertex) override {
                
                auto cluster = std::static_pointer_cast<typename Final::Graph>(g);
                auto prop = cluster->template getProperty<GeometryProperty>(lv);
                auto functor = PropRetriever(prop);
                InheritedV::apply(functor);
            };
            
        private:
            struct PropCreator : dcm::visitor<symbolic::Geometry*> {              
                template<typename T>
                result_type operator()(const T&) const {
                    typedef geometry::extractor<typename geometry_traits<T>::type> extractor;
                    typedef symbolic::TypeGeometry<typename extractor::template primitive<typename Final::Kernel>> TypeGeometry;
                    return new TypeGeometry();
                };
            };
            
            struct PropAssigner : dcm::visitor<>{  
                symbolic::Geometry* geom;                
                PropAssigner(symbolic::Geometry* g) : geom(g) {};
                
                template<typename T>
                result_type operator()(const T& geometry) const {
                    typedef typename Final::Kernel::Scalar Scalar;
                    typedef geometry::extractor<typename geometry_traits<T>::type> extractor;
                    typedef symbolic::TypeGeometry<typename extractor::template primitive<typename Final::Kernel>> TypeGeometry;

                    (typename geometry_traits<T>::modell()).template extract<Scalar, typename geometry_traits<T>::accessor >(geometry, 
                                                                                static_cast<TypeGeometry*>(geom)->getPrimitve());
                };
            };
            
            struct PropRetriever : dcm::visitor<> {  
                symbolic::Geometry* geom;                
                PropRetriever(symbolic::Geometry* g) : geom(g) {};
                
                template<typename T>
                void operator()(T& geometry) const {
                    typedef typename Final::Kernel::Scalar Scalar;
                    typedef geometry::extractor<typename geometry_traits<T>::type> extractor;
                    typedef symbolic::TypeGeometry<typename extractor::template primitive<typename Final::Kernel>> TypeGeometry;

                    (typename geometry_traits<T>::modell()).template inject<Scalar, typename geometry_traits<T>::accessor >(geometry, 
                                                                                static_cast<TypeGeometry*>(geom)->getPrimitve());
                };
            };
        };
        
        struct Constraint3D : public Stacked::Object {
            
            typedef utilities::Variant<types...> InheritedV;
            typedef typename Stacked::Object     InheritedO;
                       
            //typedef symbolic::ConstraintProperty ConstraintProperty;
            typedef symbolic::ConstraintListProperty ConstraintList;
            
            DCM_OBJECT_ADD_PROPERTIES( InheritedO, (ConstraintList) )
            
        public:
            Constraint3D(Final* system) 
                : Stacked::Object ( Final::template objectTypeID<typename Final::Geometry3D>::ID::value ),
                m_system(system){
                    
            };    
            
            template<typename ...Constraints>
            void setConstraints(Constraints&... cons) {
                
                auto cluster = std::static_pointer_cast<typename Final::Graph>(m_system->getGraph());
                
                if(m_geometries.empty())
                    throw creation_error() <<  boost::errinfo_errno(21) << error_message("Geometry was not set before setting constraints types");
                
                std::vector<symbolic::Constraint*> vec = getConstraintList();
               
                //ensure we do not have any constraints left in the graph
                if(m_geometries.size() == 1) {
                    auto& cons = cluster->template getPropertyAccessible<symbolic::ConstraintListProperty>(m_geometries[0]->getVertexProperty());
                    //remove all elemets that are in both vectors from vector 1
                    auto end = std::set_difference(cons.begin(), cons.end(), vec.begin(), vec.end(), cons.begin());
                    cons.erase(end, cons.end());
                } else {                    
                    auto edges = cluster->edges(m_geometries[0]->getVertexProperty(), m_geometries[1]->getVertexProperty());
                    for(const graph::GlobalEdge& edge : edges) {
                        symbolic::Constraint* c = cluster->template getProperty<symbolic::ConstraintProperty>(edge);
                        if(std::find(vec.begin(), vec.end(), c) != vec.end()) {
                                delete c;
                                cluster->removeEdge(edge);
                        }
                    }
                }
                                             
                //create a new global edge for every constraint
                vec.clear();
                expand(createGraphRepresentation(cons, vec)...);                
                setConstraintList(vec);
            };
            
            void set(std::shared_ptr<Geometry3D> G1, std::shared_ptr<Geometry3D> G2) {
                
                //it makes no sense to allow a geometry change, the user shall rather delete this and create a new
                //constraint
                if(!m_geometries.empty())
                    throw creation_error() <<  boost::errinfo_errno(23) << error_message("Geometries for this constraint are already set");
                
                if(!G1->holdsGeometry() || !G2->holdsGeometry())
                    throw creation_error() <<  boost::errinfo_errno(26) << error_message("Geometry must be valid");
                
                m_geometries.resize(2);
                m_geometries[0] = G1;
                m_geometries[1] = G2;                
            };
            
            //for single geometry constraint
            void set(std::shared_ptr<Geometry3D> G) {
                
                //it makes no sense to allow a geometry change, the user shall rather delete this and create a new
                //constraint
                if(!m_geometries.empty())
                    throw creation_error() <<  boost::errinfo_errno(23) << error_message("Geometry for this constraint are already set");
                
                if(!G->holdsGeometry())
                    throw creation_error() <<  boost::errinfo_errno(26) << error_message("Geometry must be valid");
                
                m_geometries.resize(1);
                m_geometries[0] = G;
            };
            
        protected:
            Final* m_system;
            std::vector<std::shared_ptr<Geometry3D>> m_geometries;
            
            template<typename ...Args>
            void expand(const Args&... args){};
            
            template<typename T>
            int createGraphRepresentation(T& t, std::vector<symbolic::Constraint*>& vec) {
                
                if(T::Arity != m_geometries.size())
                    throw creation_error() <<  boost::errinfo_errno(25) << error_message("Constraint does not support given amount of geometries");
               
                symbolic::TypeConstraint<T>* tc;
                
                //single geometry constraints are created in a special manner
                if(T::Arity == 1) {
                    //add the primitive constraint to the vertex
                    tc = new symbolic::TypeConstraint<T>();
                    tc->setPrimitive(t);
                    
                    std::shared_ptr<typename Final::Graph> cluster = std::static_pointer_cast<typename Final::Graph>(m_system->getGraph());
                    auto& cons = cluster->template getPropertyAccessible<symbolic::ConstraintListProperty>(m_geometries[0]->getVertexProperty());
                    cons.push_back(tc);
                }
                else if(T::Arity == 2) {
        
                    //let's create a new global edge for this constraint
                    std::shared_ptr<typename Final::Graph> cluster = std::static_pointer_cast<typename Final::Graph>(m_system->getGraph());
                    fusion::vector<graph::LocalEdge, graph::GlobalEdge, bool, bool> res = 
                        cluster->addEdge(m_geometries[0]->getVertexProperty(), m_geometries[1]->getVertexProperty());
                        
                    //check if we have been successfull
                    if(!fusion::at_c<2>(res))
                        throw creation_error() <<  boost::errinfo_errno(22) << error_message("Graph representation of constraint could not be created");                
                    
                    //add the primitive constraint to the global edge
                    tc = new symbolic::TypeConstraint<T>();
                    tc->setPrimitive(t);
                    cluster->template setProperty<symbolic::ConstraintProperty>(fusion::at_c<1>(res), tc);
                }
                else 
                    throw creation_error() <<  boost::errinfo_errno(25) << error_message("Constraint does not support given amount of geometries");
            
                t.setDefault();               
                vec.push_back(tc);
                
                //return a type needed for expand function to allow to call this function on parameter packs
                return 0;
            };
        };
        
        DCM_MODULE_ADD_OBJECTS(Stacked, (Geometry3D)(Constraint3D))
        
        template<typename T>
        std::shared_ptr<Geometry3D> addGeometry3D(const T& geom) {
            
            auto g = std::make_shared<typename Final::Geometry3D>(static_cast<Final*>(this));
            g->set(geom);
            return g;
        };
        
        template<typename ...Constraints>
        std::shared_ptr<Constraint3D> addConstraint3D(std::shared_ptr<Geometry3D> G1,
                                                      std::shared_ptr<Geometry3D> G2,
                                                      Constraints&... cons) { 
            
            auto c = std::make_shared<typename Final::Constraint3D>(static_cast<Final*>(this));
            c->set(G1, G2);
            c->setConstraints(cons...);
            return c;
        };
        
        template<typename ...Constraints>
        std::shared_ptr<Constraint3D> addConstraint3D(std::shared_ptr<Geometry3D> G,
                                                                      Constraints&... cons) { 
            
            auto c = std::make_shared<typename Final::Constraint3D>(static_cast<Final*>(this));
            c->set(G);
            c->setConstraints(cons...);
            return c;
        };
        
    protected:
        friend struct Geometry3D;
        friend struct Constraint3D;
        
        DCM_MODULE_ADD_GEOMETRIES(Stacked, (geometry::Point3)(geometry::Line3)(geometry::Plane)(geometry::Cylinder)
                                           (geometry::Cluster3))
    };
    
};

}//dcm

#endif //DCM_GEOMETRY3D_H







