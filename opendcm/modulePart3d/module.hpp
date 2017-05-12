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
       
enum class Scope {Local, Part, Global};
    
template<typename Type>
struct ModulePart3D {

    typedef boost::mpl::int_<6> ID;

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
            
            DCM_OBJECT_ADD_PROPERTIES( Final, (VertexProperty) )
                        
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
                cluster->template setProperty<details::GraphObjectProperty>(res.second, Inherited::shared_from_this());
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
                cluster->template setProperty<details::GraphObjectProperty>(res.second, Inherited::shared_from_this());
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
            };
            
            /**
             * @brief Adds geometry to the Part
             * The geometry is rigidly connected with all other geometries within the Part, that means
             * all constraints with other inner part geometries are ignored. However, constraints to 
             * geometries outside the Part or within other Parts can be created. They are than resolved 
             * by moving the Part in such a way, that the geometry constraints are fullfiled.
             * @note when adding geometrie it is important in which coordinate system it is defined. 
             *       This can be specified with the second parameter, there are 3 options:
             *       1. Local, meaning the geometry is defined in the Parts coordinate system
             *       2. Part, meaning the geometry is specified in the same CS as the part itself.
             *       3. Global, meaning the geometry is specified within the toplevel reference CS. If 
             *          the Part is not within annother Part "Global" is equal to "Part", but if the 
             *          the Part is a subcomponent of other Parts all of them will be taken into account.
             */
            template<typename T>
            std::shared_ptr<typename Stacked::Geometry3D> addGeometry3D(const T& geom, Scope s = Scope::Local) {
                
                auto g = m_system->addGeometry3D(geom);
                
                //we need to transfer the geometry to our cluster
                auto cluster = std::static_pointer_cast<typename Final::Graph>(m_system->getGraph());
                auto globalVertex = g->template getProperty<VertexProperty>();
                //it may be a sub-sub cluster etc.
                auto localVertex = cluster->getLocalVertex(globalVertex);
                dcm_assert(localVertex.second);
                cluster->moveToSubcluster(localVertex.first, getVertexProperty());    
                
                //store the scope for later transformations
                m_scopeMap[g] = s;
                
                return g;
            };
            
            /**
             * @brief Adds a child Part to this Part
             * The new part is rigidly connected with all other geometries within the Part, that means
             * all constraints with other inner part geometries are ignored. However, constraints to 
             * geometries outside the Part or within other Parts can be created. They are than resolved 
             * by moving the parent Part in such a way, that the child Part constraints are fullfiled.
             * @note when adding a part it is important in which coordinate system it is defined. 
             *       This can be specified with the second parameter, there are 3 options:
             *       1. Local, meaning the part is defined in the Parts coordinate system
             *       2. Part, meaning the part is specified in the same CS as the part itself.
             *       3. Global, meaning the part is specified within the toplevel reference CS. If 
             *          the parent Part is not within annother Part "Global" is equal to "Part", but if the 
             *          the parent Part is a subcomponent of other Parts all of them will be taken into account.
             */
            std::shared_ptr<typename Stacked::Geometry3D> addPart3D(const Type& part, Scope s = Scope::Local) {
                                
                auto p = std::make_shared<Part3D>(std::static_pointer_cast<Part3D>(Inherited::shared_from_this()));
                p->set(part);
                
                //store the scope for later transformations
                m_scopeMap[p] = s;
                
                return p;
            };
            
            bool holdsGeometry() {return m_isSet;}
            
        protected:
            Final*                                                      m_system;
            Type                                                        m_type;
            bool                                                        m_isSet = false;
            std::shared_ptr<graph::AccessGraphBase>                     m_cluster;
            std::map<std::shared_ptr<dcm::details::GraphObject>, Scope> m_scopeMap;
            
            virtual void preprocessVertex(std::shared_ptr<graph::AccessGraphBase> g, 
                                          graph::LocalVertex lv, graph::GlobalVertex gv) override {
                
                auto cluster = std::static_pointer_cast<typename Final::Graph>(g);
                auto prop = cluster->template getProperty<GeometryProperty>(lv);
                
                //if the property is not available or does not hold the correct type we need to delete it, no way around
                if(prop && prop->getType() != geometry::Part3<Kernel>::index()) {
                    delete prop;
                    prop = nullptr;
                }
                
                //if no property in existance we need to create it
                if(!prop) {
                    prop = new symbolic::TypeGeometry<geometry::Part3<Kernel>>;
                    cluster->template setProperty<GeometryProperty>(lv, prop);
                }
                
                //TODO: transform the child geometries dependend on their scope
                    
                //we definitly need to set the value
                (typename geometry_traits<Type>::modell()).template extract<typename Kernel::Scalar, 
                                                                    typename geometry_traits<Type>::accessor>(
                              m_type, static_cast<symbolic::TypeGeometry<geometry::Part3<Kernel>>*>(prop)->getPrimitve());
            };
            
            virtual void postprocessVertex(std::shared_ptr<graph::AccessGraphBase> g,
                                           graph::LocalVertex lv, graph::GlobalVertex) override {
                
                auto cluster = std::static_pointer_cast<typename Final::Graph>(g);
                auto prop = cluster->template getProperty<GeometryProperty>(lv);
                (typename geometry_traits<Type>::modell()).template inject<typename Kernel::Scalar, 
                                                                    typename geometry_traits<Type>::accessor >(
                             m_type,  static_cast<symbolic::TypeGeometry<geometry::Part3<Kernel>>*>(prop)->getPrimitve());
            };
        };
        
        DCM_MODULE_ADD_OBJECTS(Stacked, (Part3D))

    public:  
        std::shared_ptr<Part3D> addPart3D(const Type& part) {
            
            auto p = std::make_shared<Part3D>(static_cast<Final*>(this));
            p->set(part);
            return p;
        };
        
    protected:
        friend struct Part3D;
        
    };
    
};

}//dcm

#endif //DCM_GEOMETRY3D_H







