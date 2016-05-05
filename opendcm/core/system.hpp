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

#ifndef DCM_SYSTEM_H
#define DCM_SYSTEM_H

#include <boost/mpl/less_equal.hpp>
#include <boost/mpl/fold.hpp>
#include <boost/mpl/sort.hpp>
#include <boost/mpl/vector.hpp>

#include <boost/preprocessor/enum_params.hpp>
#include <boost/preprocessor/enum_params_with_a_default.hpp>

#include "module.hpp"
#include "kernel.hpp"

#define DCM_MAX_MODULE_SIZE 5

namespace mpl = boost::mpl;

namespace dcm {
namespace details {

template<typename M1, typename M2>
struct module_sort : mpl::less_equal<typename M1::ID, typename M2::ID> {};

template<typename Final, typename Stack, typename Module>
struct module_inheritance {
    typedef typename Module::template type<Final, Stack> type;
};

template<typename Final, typename ...ModuleTypes>
struct module_inheriter {

    typedef mpl::vector<ModuleTypes...> Modules;
    
    //find or create the kernel
    typedef typename mpl::find_if<Modules, boost::is_base_of<numeric::KernelBase, mpl::_1> >::type it;
    typedef typename mpl::if_<boost::is_same<it, typename mpl::end<Modules>::type>, Eigen3Kernel<double>, 
                                typename mpl::deref<it>::type>::type Kernel;   
    
    //remove kernel from modules list if provided
    typedef typename mpl::if_<boost::is_same<it, typename mpl::end<Modules>::type>, Modules, 
                                typename mpl::remove<Modules,Kernel>::type>::type StrippedModules;    
    
    //sort according to the modules id
    typedef typename mpl::sort<StrippedModules, module_sort<mpl::_1, mpl::_2> >::type SortedModules;

    //initialise the module stack
    typedef typename mpl::fold<SortedModules, ModuleCoreInit<Final, Kernel>,
            module_inheritance<Final, mpl::_1, mpl::_2> >::type ModuleStack;

    //and create the finished inheritance type
    typedef ModuleCoreFinish<Final, ModuleStack > type;
};


}; //details

template<typename ...Modules>
class System : public details::module_inheriter< System<Modules...>, Modules... >::type {};

}; //dcm

#endif //GCM_SYSTEM_H













