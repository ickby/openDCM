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
    GNU Lesser General Public License for more detemplate tails.

    You should have received a copy of the GNU Lesser General Public License along
    with this library; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/

#ifndef DCM_CONSTRAINT_H
#define DCM_CONSTRAINT_H

#include <boost/utility/enable_if.hpp>
#include <boost/mpl/contains.hpp>
#include <boost/fusion/include/vector.hpp>
#include <boost/fusion/include/find.hpp>
#include <boost/fusion/include/at.hpp>

#include <iostream>

namespace mpl = boost::mpl;
namespace fusion = boost::fusion;

namespace dcm { 
    
//the possible directions
enum class Orientations { Parallel, Equal, Opposite, Perpendicular };

//the possible solution spaces
enum class SolutionSpaces {Bidirectional, Positiv_directional, Negative_directional};

/**
 * @brief Constraint primitives handling
 *
 * Equivalent to the geometric primitives there is a definition for constraint primitives. A constraint 
 * primitive is a constraint type like distace, primitive in this context means that it holds only its type
 * and all needed information to fully define it, called option. For a distance constraint this is for 
 * example the nuemeric distance value. A primitive is used to store and access the type and its data in a
 * compfortable and user friendly way.
 * There are a few requirements for primitive constraints. First they must be assignalble and copyconstructible
 * while preserving their option values. Second all options need to be stored in a  fusion vector called
 * m_storage. Third it must be possible to provide options via operator() for one or many options at once 
 * and furthermore through assigning with operator=, also for single options or multiple ones via initializer
 * lists. The last requirement regards the default values of the options. It is important to have the possibility
 * to reset all options to default during the objects lifetime. Therefore the function setDefault() must exist.
 * An example primitive constraint which holds all requirements could look like this:
 * @code
 * struct TestConstraint {
 *   
 *  fucion::vector<int, char> m_sequence; //requirement 2
 *  
 *  //requirement 3
 *  void operator()(int);
 *  void operator()(char);
 *  void operator()(int, char);
 *  void operator=(int);
 *  void operator=(char);
 *  void operator=(int, char);
 * 
 *  void setDefault(); //requirement 4;
 * };
 * @endcode
 * 
 * Requirement 1 is fullfilled by the standart compiler generated assignement operator and copy constructor. This
 * also holds for the assignemnd of multiple options via initializer lists.
 * 
 * \Note The multiassignement via operator() must not support arbitrary argument ordering. It is sufficient to 
 * allow the parameter only in the same order as the options are specified.
 */
namespace constraint {  

    /**
     * @brief Herlper class to create primitive constraints
     * 
     * The primitive constraint type concept formulates a few requirements on the used types. To achieve a full
     * compatibility with the concept a certain boilerplate is needed. To circumvent this boilerplate this class
     * is given. It accepts the option types for the constraint as template arguments and then provides
     * all needed storages and assignment operators. 
     * This only works for option types which are default constructible. The only second restriction is that the
     * option types need to distuinguishable. To use a option type twice is not allowd, for example two times int.
     * Otherwise the option assignement of single options would be ambigious. 
     * The default value requirement is achieved through the constructor which acceppts the default values. If 
     * the standart constructor is used then the default constructed values of the option types are used as 
     * default values. Note that you must provide the constructors also in your derived class.
     * An example cosntraint can look like this:
     * @code
     * struct TetstConstraint : public constraint::Constraint<int, char> {
     *    TestConstraint(const int& i, const char& c) : Constraint(i,c) {};
     * };
     * TestConstraint test(1,'a');
     * @endcode
     * 
     * \param OptionTypes a variadic sequence of copy constructable types which describe the stored options
     */    
template<typename ...OptionTypes>
struct Constraint {

    typedef typename fusion::vector<OptionTypes...> Options;

    Constraint() {};
    
    Constraint(const OptionTypes&... defaults) : m_defaults(defaults...) {
        m_storage = m_defaults;
    };
    
    //Copy assign option. We provide this to allow the automatic copy constructor to be generated. 
    //This gives 2 operators for the price of implementing one, for this base and all derived classes.
    //As we have only one parameter in this function it would be ambigious with the single option 
    //operator(), hence we need to activate/deactivate both according to the provided type.
    template<typename T>
    typename boost::disable_if<mpl::contains<Options, T>, Constraint&>::type operator()(const T& c) {
        m_storage = c.m_storage;
        return *this;
    };

    //Set a single option. As we have only one parameter in this function it would be ambigious with the 
    //copy operator(), hence we need to activate/deactivate both according to the provided type.
    template<typename T>
    typename boost::enable_if<mpl::contains<Options, T>, Constraint&>::type operator()(const T& val) {
        BOOST_MPL_ASSERT((mpl::contains<Options, T>));
        *fusion::find<T>(m_storage) = val;
        return *this;
    };
    
    //set multiple options at once.
    Constraint& operator()(const OptionTypes&... val) {
        m_storage = Options(val...);
        return *this;
    };

    //Assign option. Disable it to avoid confusion with the copy assignement operator
    template<typename T>
    typename boost::enable_if<mpl::contains<Options, T>, Constraint&>::type operator=(const T& val) {
        return operator()(val);
    };
    
    //explicity give copy assignement as it otherwise would be treated as deleted due to the assignend operator
    //given above. When this one is provided the derived class wil also generate one.
    Constraint& operator=(const Constraint& val) {
        m_storage = val.m_storage;
        return *this;
    };
   
    //set default option values, neeeded for repedability and to prevent unexpected behaviour
    void setDefault() {
        m_storage = m_defaults;
    };
       
protected:
    Options       m_storage;
    const Options m_defaults;
};

}//constraint
    
namespace numeric {
    

/**
 * @brief Basic numeric class for evaluating primitive constraints
 * 
 * Primitive constraints hold only their constraint type and the options to fully define its behaviour. This is
 * however not enough to calculate the error functions and jacobi entries in a numerical solving context as it 
 * does not provide any equations. 
 * 
 * \param Kernel the math kernel in use
 * \param PC the primitive constraint in use
 */
template<typename Kernel, typename PC, template<class, bool> class PG1, template<class, bool> class PG2>
struct Constraint {
    
private:
    //numeric::Geometry G1<Kernel, PG1> m_geometry1;
    //numeric::Geometry G2<Kernel, PG2> m_geometry2;
};
    
}//numeric
    
    
namespace symbolic {
    
struct Constraint {
    
    int type;
};

template<typename C>
struct TypeConstraint : public Constraint {

    typedef C PrimitiveConstraint; 
    
    PrimitiveConstraint& getPrimitveConstraint() {
        return m_constraint;
    };
    
    void setPrimitiveConstraint(const C& c) {
        //type = id;
        m_constraint = c;
    };
    
    void setConstraintID(int id) {
        type = id;
    }
    
protected:   
    PrimitiveConstraint m_constraint;
};

struct ConstraintProperty {
    typedef Constraint* type;
    struct default_value {
        Constraint* operator()() {
            return nullptr;
        };
    };
    struct change_tracking{};
};
    
}//symbolic    



//Provide a few default constraints.
//static is needed to restrain the scope of the objects to the current compilation unit. Without it
//every compiled file including this header would define these as global and the linker would find
//multiple definitions of the same objects

struct Distance : public dcm::constraint::Constraint<double, SolutionSpaces> {
    using Constraint::operator=;
    Distance(const double& i, SolutionSpaces s) : Constraint(i,s) {};
    
    double&        value() {return fusion::at_c<0>(m_storage);};
    SolutionSpaces solutionSpace() {return fusion::at_c<1>(m_storage);};
};

struct Orientation : public dcm::constraint::Constraint<Orientations> {
    using Constraint::operator=;
    Orientation(const Orientations& i) : Constraint(i) {};
    
    Orientations& value() {return fusion::at_c<0>(m_storage);};
};

struct Angle : public dcm::constraint::Constraint<double> {
    using Constraint::operator=;
    Angle(const double& i) : Constraint(i) {};
    
    double& value() {return fusion::at_c<0>(m_storage);};
};

static Distance         distance(0, SolutionSpaces::Bidirectional);
static Orientation      orientation(Orientations::Parallel);
static Angle            angle(0);

}//dcm

#endif //DCM_CONSTRAINT_H





