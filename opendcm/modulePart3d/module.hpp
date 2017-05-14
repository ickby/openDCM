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

#ifndef DCM_MODULE_PART3D_H
#define DCM_MODULE_PART3D_H

#include "opendcm/core/object.hpp"
#include "opendcm/core/module.hpp"
#include "opendcm/core/utilities.hpp"
#include "opendcm/core/geometry.hpp"
#include "opendcm/core/constraint.hpp"
#include "opendcm/core/typeadaption.hpp"

#include "geometry.hpp"


#include <type_traits>

namespace dcm {
       
enum class Scope {Local, Global};
    
template<typename Type>
struct ModulePart3D {

    typedef boost::mpl::int_<6> ID;

    template<typename Final, typename Stacked>
    struct type : public Stacked {

        typedef typename Stacked::Kernel Kernel;
        typedef dcm::details::Transform<typename Kernel::Scalar, 3> Transform;
        
        //setup a transform property for stacked transforms
        struct Transform3Property {
            typedef Transform type;
        };
        
        DCM_MODULE_ADD_CLUSTER_PROPERTIES(Stacked, (Transform3Property));
        
        
        type() : Stacked() {};        
        ~type() {}
        
        //we setup the reduction graph after the constructors are done, as only then the graph 
        //is available (as well as the virtual function to get it)
        void init() {
            
            //that is our responsibility in the init stack
            Stacked::init();
            
        };
        
        struct Part3D;
        
        struct Geometry3D : public Stacked::Geometry3D {
            
            Geometry3D(Final* system) : Stacked::Geometry3D(system) {};
                
            std::shared_ptr<Part3D> parentPart() {
                if(!m_parent.expired())
                    return m_parent.lock();
                return std::shared_ptr<Part3D>();
            };
            
        protected:
            typedef typename Stacked::Geometry3D Base;

            std::weak_ptr<Part3D> m_parent;
                        
            void transform(const Transform& trans) {
            
                //to access the right primitive type and transform it we use the visitors
                symbolic::Geometry* prop;
                if(!m_parent.expired()) 
                    prop = std::static_pointer_cast<typename Final::Graph>(m_parent.lock()->m_cluster)->template getProperty<symbolic::GeometryProperty>(Base::getVertexProperty());
                else {
                    auto cluster = std::static_pointer_cast<typename Final::Graph>(Base::m_system->getGraph());
                    prop = cluster->template getProperty<symbolic::GeometryProperty>(Base::getVertexProperty());
                }
                PropTransformer functor(prop, trans);
                Base::apply(functor);
            };
            
            friend struct Part3D;
            
        private:
            struct PropTransformer : dcm::visitor<> {  
                symbolic::Geometry* m_geom;            
                const Transform&    m_trans;
                
                PropTransformer(symbolic::Geometry* g,
                              const Transform& t) : m_geom(g), m_trans(t) {};
                
                template<typename T>
                void operator()(T& geometry) const {
                    typedef typename Final::Kernel::Scalar Scalar;
                    typedef geometry::extractor<typename geometry_traits<T>::type> extractor;
                    typedef symbolic::TypeGeometry<typename extractor::template primitive<typename Final::Kernel>> TypeGeometry;

                    static_cast<TypeGeometry*>(m_geom)->getPrimitve().transform(m_trans);
                };
            };
        };
        
        
        /**
         * @brief Container for 3D user parts
         * 
         * A Part is a container for geometry, which is rigidly connected. It provides a way to move 
         * the rigid system in 6 degrees of freedom. This movement is what defines the Part and what 
         * the user must supply via the given Module type.
         * 
         * To add geometry to the part simple use the addGeometry() method of the object excactly the 
         * same as its system counterpart. 
         */
        struct Part3D : public Stacked::Object {

            typedef typename Stacked::Object   Inherited;   
            
            typedef symbolic::GeometryProperty GeometryProperty;
            typedef graph::VertexProperty      VertexProperty;
            
            DCM_OBJECT_ADD_PROPERTIES( Inherited, (VertexProperty) )
                        
        public:            
            Part3D(Final* system) 
                : Stacked::Object( Final::template objectTypeID<typename Final::Part3D>::ID::value ),
                m_system(system) {
                    
                typedef typename Final::Kernel  Kernel;
                typedef typename Kernel::Scalar Scalar;
                typedef symbolic::TypeGeometry<Part3D> TypeGeometry;
                
                //we may need to setup the graph. We do this here to allow to clear the geometry and later reinitialize by setting a new geometry.
                std::shared_ptr<typename Final::Graph> cluster = std::static_pointer_cast<typename Final::Graph>(m_system->getGraph());
                std::pair<std::shared_ptr<typename Final::Graph>, graph::LocalVertex> res = cluster->createCluster();
                setVertexProperty(cluster->getGlobalVertex(res.second));
                m_cluster = res.first;
            
            };
            
            Part3D(std::shared_ptr<Part3D> parent) 
                : Stacked::Object( Final::template objectTypeID<typename Final::Part3D>::ID::value ),
                m_system(parent->m_system) {
                    
                typedef typename Final::Kernel  Kernel;
                typedef typename Kernel::Scalar Scalar;
                typedef symbolic::TypeGeometry<Part3D> TypeGeometry;
                
                //we may need to setup the graph. We do this here to allow to clear the geometry and later reinitialize by setting a new geometry.
                std::shared_ptr<typename Final::Graph> cluster = std::static_pointer_cast<typename Final::Graph>(parent->m_cluster);
                std::pair<std::shared_ptr<typename Final::Graph>, graph::LocalVertex> res = cluster->createCluster();
                setVertexProperty(cluster->getGlobalVertex(res.second));
                m_cluster = res.first;
            };
            
            /**
             * @brief Set the content and stores the user type
             * 
             * This function is used to set the content of the Part3D container. It extract the contents of 
             * the given user Type and stores it in an accessible manner for the solver. Furhtermore it stores
             * the supplied user type for later access. The exact behaviour of this feature depends on the specified
             * qualifiers of the storable types as given as the Module3D template arguments. If the type was given
             * without any special qualifiers then a copy of the given user type is strored.
            */
            void set(const Type& part) {
                //store the type
                m_type = part;
                m_isSet  = true;
                
                // we were not able to set the GraphObject before, as shared_from_this was not available in the constructor
                std::shared_ptr<typename Final::Graph> cluster = std::static_pointer_cast<typename Final::Graph>(m_system->getGraph());
                auto res = cluster->getLocalVertexGraph(getVertexProperty()); //fusion::vector<LocalVertex, std::shared_ptr<ClusterGraph>, bool>
                dcm_assert(fusion::at_c<2>(res));
                fusion::at_c<1>(res)->template setProperty<details::GraphObjectProperty>(fusion::at_c<0>(res), Inherited::shared_from_this());
            };
            
            /**
             * @brief Get the content as usertype
            */
            const Type& get() {
                //store the type
                return m_type;
            };
            
            /**
             * @brief Adds geometry to the Part
             * The geometry is rigidly connected with all other geometries within the Part, that means
             * all constraints with other inner part geometries are ignored. However, constraints to 
             * geometries outside the Part or within other Parts can be created. They are than resolved 
             * by moving the Part in such a way, that the geometry constraints are fullfiled.
             * @note when adding geometrie it is important in which coordinate system it is defined. 
             *       This can be specified with the second parameter, there are these options:
             *       1. Local, meaning the geometry is defined in the Parts coordinate system
             *       2. Global, meaning the geometry is specified within the toplevel reference CS. If 
             *          the Part is not within annother Part "Global" is equal to "Part", but if the 
             *          the Part is a subcomponent of other Parts all of them will be taken into account.
             */
            template<typename T>
            std::shared_ptr<Geometry3D> addGeometry3D(const T& geom, Scope s = Scope::Local) {
                
                auto g = m_system->addGeometry3D(geom);
                
                //we need to transfer the geometry to our cluster
                auto cluster = std::static_pointer_cast<typename Final::Graph>(m_system->getGraph());
                auto globalVertex = g->getVertexProperty();
                //it may be a sub-sub cluster etc.
                auto localVertex = cluster->getLocalVertex(globalVertex);
                dcm_assert(localVertex.second);
                cluster->moveToSubcluster(localVertex.first, getVertexProperty());    
                
                //store the scope and geometry for later transformations
                m_geometries.push_back(std::make_pair(g, s));
                
                return g;
            };
            
            /**
             * @brief Adds a child Part to this Part
             * The new part is rigidly connected with all other geometries within the Part, that means
             * all constraints with other inner part geometries are ignored. However, constraints to 
             * geometries outside the Part or within other Parts can be created. They are than resolved 
             * by moving the parent Part in such a way, that the child Part constraints are fullfiled.
             * @note when adding a part it is important in which coordinate system it is defined. 
             *       This can be specified with the second parameter, there are these options:
             *       1. Local, meaning the part is defined in the Parts coordinate system
             *       2. Global, meaning the part is specified within the toplevel reference CS. If 
             *          the parent Part is not within annother Part "Global" is equal to "Part", but if the 
             *          the parent Part is a subcomponent of other Parts all of them will be taken into account.
             */
            std::shared_ptr<typename Stacked::Geometry3D> addPart3D(const Type& part, Scope s = Scope::Local) {
                                
                auto p = std::make_shared<Part3D>(std::static_pointer_cast<Part3D>(Inherited::shared_from_this()));
                p->set(part);
                
                //store the scope for later transformations
                p->m_scope = s;
                
                return p;
            };
            
            bool holdsGeometry() {return m_isSet;}
            
        protected:
            Final*                                                      m_system;
            Type                                                        m_type;
            bool                                                        m_isSet = false;
            Scope                                                       m_scope;
            Transform                                                   m_transform;
            std::shared_ptr<graph::AccessGraphBase>                     m_cluster;
            std::vector<std::pair<std::shared_ptr<Geometry3D>,Scope>>   m_geometries;
            
            virtual void preprocessCluster(std::shared_ptr<graph::AccessGraphBase> g, graph::LocalVertex v, 
                                           std::shared_ptr<graph::AccessGraphBase> lg) override {
            
                //Set the value for the Part. We need to do it here, and not in preprocessVertex, 
                //as we need the transform to calculate the stacked transform
                auto cluster = std::static_pointer_cast<typename Final::Graph>(g);
                auto prop = cluster->template getProperty<GeometryProperty>(v);
                
                //if the property is not available or does not hold the correct type we need to delete it, no way around
                if(prop && prop->getType() != geometry::Part3<Kernel>::id()) {
                    delete prop;
                    prop = nullptr;
                }
                
                //if no property in existance we need to create it
                if(!prop) {
                    prop = new symbolic::TypeGeometry<geometry::Part3<Kernel>>;
                    cluster->template setProperty<GeometryProperty>(v, prop);
                }
                                   
                //we definitly need to set the value
                auto& primitive = static_cast<symbolic::TypeGeometry<geometry::Part3<Kernel>>*>(prop)->getPrimitve();            
                
                (typename geometry_traits<Type>::modell()).template extract<typename Kernel::Scalar, 
                                                                    typename geometry_traits<Type>::accessor>(
                                                                                m_type, primitive);
                                                                    
                //in dcm everything must be defined in local scope!
                if(m_scope == Scope::Global)
                    primitive.transform(primitive.transform().inverse()*cluster->template getProperty<Transform3Property>());
      
                
                //calculate the global transforms
                m_transform = primitive.transform();
                std::static_pointer_cast<typename Final::Graph>(lg)->template setProperty<Transform3Property>(
                     cluster->template getProperty<Transform3Property>()*primitive.transform());   
            }
            
            virtual void preprocessVertex(std::shared_ptr<graph::AccessGraphBase> g, 
                                          graph::LocalVertex lv, graph::GlobalVertex gv) override {
                
                //now all vertices have been initialised. We can transform them to be global! 
                for(auto gpair : m_geometries) {                    
                    if(gpair.second == Scope::Global) {
                        auto& trans = std::static_pointer_cast<typename Final::Graph>(m_cluster)->template getProperty<Transform3Property>();
                        gpair.first->transform(trans.inverse());
                    }
                }
            };
            
            virtual void postprocessCluster(std::shared_ptr<graph::AccessGraphBase> g, graph::LocalVertex v, 
                                            std::shared_ptr<graph::AccessGraphBase> lg) override {
              
                auto cluster = std::static_pointer_cast<typename Final::Graph>(g);
                auto prop = cluster->template getProperty<GeometryProperty>(v);
                auto& primitive = static_cast<symbolic::TypeGeometry<geometry::Part3<Kernel>>*>(prop)->getPrimitve();
                
                //set the value
                (typename geometry_traits<Type>::modell()).template inject<typename Kernel::Scalar, 
                                                                    typename geometry_traits<Type>::accessor >(
                                                                            m_type,  primitive);
                                                                                    
                //calculate the transforms
                m_transform = primitive.transform();
                std::static_pointer_cast<typename Final::Graph>(lg)->template setProperty<Transform3Property>(
                            cluster->template getProperty<Transform3Property>()*primitive.transform());
                
                //set it as the user expects it
                if(m_scope == Scope::Global)
                    primitive.transform(std::static_pointer_cast<typename Final::Graph>(g)->template getProperty<Transform3Property>());
            };
                
            virtual void postprocessVertex(std::shared_ptr<graph::AccessGraphBase> g,
                                           graph::LocalVertex lv, graph::GlobalVertex) override {
                
                //now all vertices have been initialised. We can transform them to be global! 
                for(auto gpair : m_geometries) {                    
                    if(gpair.second == Scope::Global) {
                        auto& trans = std::static_pointer_cast<typename Final::Graph>(m_cluster)->template getProperty<Transform3Property>();
                        gpair.first->transform(trans);
                    }
                }
            };
            
            friend struct Geometry3D;
        };
        
        DCM_MODULE_ADD_OBJECTS(Stacked, (Part3D))

    public:  
        std::shared_ptr<Part3D> addPart3D(const Type& part) {
            
            auto p = std::make_shared<Part3D>(static_cast<Final*>(this));
            p->set(part);
            return p;
        };
        
        //as we have created an extended Geometry3D we must override the creation function to 
        //return the latest Geometry3D type
        template<typename T>
        std::shared_ptr<Geometry3D> addGeometry3D(const T& geom) {            
            return std::static_pointer_cast<Geometry3D>(Stacked::addGeometry3D(geom));
        };
        
        //as we have redefined Geometry3D the specialisation with Geometry3D is not found by the compiler, 
        //it would always use the second version
        template<typename ...Constraints>
        std::shared_ptr<typename Stacked::Constraint3D> addConstraint3D(std::shared_ptr<Geometry3D> G,
                                                                      Constraints&... cons) { 
            
            auto c = std::make_shared<typename Final::Constraint3D>(static_cast<Final*>(this));
            c->set(G);
            c->setConstraints(cons...);
            return c;
        };
        //as we have redefined Geometry3D the specialisation with Geometry3D is not found by the compiler, 
        //it would always use this version
        template<typename ...Constraints>
        std::shared_ptr<typename Stacked::Constraint3D> addConstraint3D(std::shared_ptr<Geometry3D> G1,
                                                      std::shared_ptr<Geometry3D> G2,
                                                      Constraints&... cons) { 
            
            auto c = std::make_shared<typename Final::Constraint3D>(static_cast<Final*>(this));
            c->set(G1, G2);
            c->setConstraints(cons...);
            return c;
        };
        
        
    protected:
        friend struct Part3D;
    };
    
};

}//dcm

#endif //DCM_GEOMETRY3D_H







