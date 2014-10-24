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

/** @addtogroup Core
 * @{*/

/** @defgroup Schedule Dependency based execution framework
 *
 * @brief A possibly multithreaded framework for executiong jobs in a defined order with gurantees on 
 * when each job starts or end in relation to others. 
 * 
 * \section graph Dependencies
 * 
 * The main building block for a dependency chain is the \ref Node struct. It can hold exact one of any type 
 * of callable object which shall be executed by this \ref Node. Therefore this struct is the  equivalent of
 * a self-contained job. Inside a node no ordering can be done, when it is executed the assigned callable is
 * executed. 
 * 
 * To achieve ordering with nodes they can be dependend on each other. For example a \ref Node A can be a 
 * dependency of \ref Node B. Than it will be ensured that A is always executed before B, like A --> B. This 
 * does not mean that A must start to be executed before B but that A must Finish before B starts. Therefore
 * A can safely use all results B creates. To see how such a dependency is created have a look at the 
 * \ref Node struct.
 * 
 * Often such a simple dependency chain is not flexible enough to describe the real dependency problem. 
 * For example it can happen that an \ref Node C depends also on A, but not on B. A simple chain in the form
 * of A --> B --> C however does not describes this fact, here C will wait until B finishes even if it is 
 * unneseccary. Therefore a \ref Node can be a dependency of many other nodes. One can make C another 
 * sibling of A, something like  C <-- A ---> B. It is then guranteed that A is executed first and only after 
 * it is finished either B or C or possibly both in parallel are executed. This of course means B and C must
 * be completly independend from each other. In a multithreaded environment this also means they shall not 
 * access the same resources or make sure those resources are thread save.
 * 
 * This parallel execution of nodes is very handy, especially for multithreading. It is therefore encouraged
 * to make as much nodes parallel as possible. However, once a dependency is splitted our nodes B and C do 
 * not know when the other one is finished. So if you have a \ref Node D which depends on B and C you can 
 * simply add it as sibling of B, a chain like this C <-- A ---> B --> D would not ensure that B and C are
 * executed before D. Therefore it is possible to have multiple dependencys per node. D can simply depend on
 * B and C, something like this: C <-- A ---> B --> D <-- C. In this case it is not clear if B or C are 
 * executed first, however, it is ensured that D will only start if both are finished.
 * 
 * \section Groups
 * 
 * The described handling by nodes is very powerful. However, it also is cumbersome as you neeed to collect
 * all dependencies manual. Often one knows that a certain range of jobs that can be grouped, meaning that 
 * those jobs have any kind of dependency on each other but also need all to be done before something else
 * can happen. For this common scenario the \ref Group struct is introduces. A \ref Group is basicly a 
 * \ref Node and can be handled as one for creating dependencies. However, a \ref Group can additionaly 
 * create its own nodes. All nodes inside of a \ref Group are ensured to be started only after all of the
 * groups dependencies are finsihed. Furthermore it is guranteed that all jobs inside a \ref Group are done
 * before the groups dependencis are executed. 
 * 
 * As noted such a behaviour can be created with nodes only. However, the \ref Group struct has the big
 * advantage of enforcing the gurantees. If a \ref Node is created as part of a group one does not need to
 * ensure the gurantees manual, this is already done. Therefore it can not happen that one forgets about it.
 * 
 * As \ref Group is a \ref Node it is possible to have groups inside groups and build a certain kind if 
 * dependency-stack. 
 * 
 * \section Execution
 * 
 *@}
 */

