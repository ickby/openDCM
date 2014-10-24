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

#ifndef DCM_SHEDULER_HPP
#define DCM_SHEDULER_HPP

#include <queue>
#include <valarray>

#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <boost/exception/errinfo_errno.hpp>
#include <boost/foreach.hpp>

#ifdef DCM_USE_MULTITHREADING
#include <boost/thread.hpp>
#endif

#include "defines.hpp"

namespace dcm {
namespace details {
namespace shedule {

/**
 * @brief Exeption for sheduler errors
 *
 * This exception is thrown when the sheduler get out of sync or encounters other
 * strange errors which make a further reliable operation impossible
 **/
struct shedule_error : virtual boost::exception { };

struct Group;

/** @addtogroup Core
 * @{
 * */

/** @addtogroup Schedule
 * @{*/

/**
 * @brief Basic building block of a dependency based execution
 * 
 * 
 */
struct Node {

    Node(const boost::function<void()>& callable = NULL) 
        : registeredInput(0), handledInputs(0), m_function(callable) {};

        
    typedef std::vector<Node*>::iterator iterator;
    typedef std::vector<Node*>::const_iterator const_iterator;

    //execution path handling
    int registeredInput;
#ifdef DCM_USE_MULTITHREADING
    boost::atomic<int> handledInputs;
#else
    int handledInputs;
#endif

    /**
     * @brief Iterator to the first node that depend on \a this
     *
     * Use this to iterate over all dependend note in a standart library conformant manner
     *
     * @return NodeIterator Random access iterator pointing to the first dependend node
     * @sa endDependendNodes
     */
    inline iterator begin() {
        //the name begin() should not be changed as we then loose compatibility with boost::range 
        //and therefore with BOOST_FOREACH
        return m_next.begin();
    };

    /**
     * @brief Iterator to the past-the-end element of all dependend nodes of \a this
     *
     * Use this to iterate over all dependend note in a standart library conformant manner
     *
     * @return NodeIterator Random access iterator referring to the past-the-end element
     * @sa beginDependendNodes
     */
    inline iterator end() {
        //the name end() should not be changed as we then loose compatibility with boost::range 
        //and therefore with BOOST_FOREACH
        return m_next.end();
    };

  
    /**
     * @brief Creates a node from given callable which depends on \a this
     *
     * If a certain callable shall only be executed in case that \a this is finished the dependency can
     * be created with this function. The returned node will only execute if \a this is done. 
     * 
     * \note Even if \a this is in a group the dependend node will not ne in a group. There is no 
     * guarante when this noe will be executed in regard to any group. As long as it has no other node 
     * dependend on itself the only gurantee is that it is executed before the shedulers join() method 
     * returns
     *
     * @param node The node which depends on \a this
     * @return void
     */
    Node* createDependendNode(const boost::function<void()>& callable) {

        Node* n = new Node(callable);
        setupDependendNode(n);
        return n;
    };
    
    inline void execute() {
        if(m_function)
            m_function();
    };

protected:
    
    /**
     * @brief Connects the given node as dependend on \a this
     *
     * If a certain node shall only be executed in case that \a this is finished the dependency can
     * be setup with this function. The dependend node will only execute if \a this is done, however,
     * it must not be the directly following node that gets executed. There are no restirictions on
     * the dependend node other than it is not allowed to build up cyclic dependencies
     *
     * @param node The node which depends on \a this
     * @return void
     */
    void setupDependendNode(Node* node) {

#ifdef DCM_DEBUG
        //we don't need to impose any restrictions. A node can depend ony many others therefore
        //every possible combination appart from a circular connection is allowed. Circulars
        //however are rather costly to detect, therefore this happens for debug mode only
        dcm_assert(!node->isCyclic(this,-1));
#endif
        m_next.push_back(node);
        node->registeredInput++;
    };
    
#ifdef DCM_DEBUG
    bool isCyclic(Node* start, int depth) {
        if(depth == 0)
            return false;

        bool value = false;
        BOOST_FOREACH(Node* n, m_next) {
            value = (n==start) ? true : n->isCyclic(start, depth-1);
            if(value) 
                return true;
        }
        return false;
    };
#endif

    std::vector<Node*> m_next;
    boost::function<void()> m_function;

    friend struct Group;
};

struct Group : public Node {

    Group(const boost::function<void()>& atStart = NULL,
          const boost::function<void()>& atFinish = NULL) : Node(atStart) {

        m_final = new Node(atFinish);
    };

    ~Group() {
        delete m_final;
    };

    /**
     * @brief Creates a node in this group of connected nodes
     *
     * This function adds the callable to the execution group and ensures that it will be executed 
     * after the group setup and before the group finishes as a whole. 
     *
     * \param callable The callable object that shall be added to this execution group
     * @return The node which describes the callable inside the group
     */
    Node* createNode(const boost::function<void()>& callable) {

        Node* node = new Node(callable);

        setupGroupNode(node);        
        return node;
    };

    /**
     * @brief Creates a node in this group of connected nodes as dependend on another node
     *
     * This function adds the callable to the execution group and ensures that it will be executed after 
     * the group setup and before the group finishes as a whole. It furthermore ensures that it is executed 
     * after the finishing of \a dependsOn. It does however not ensured that \a node is executed directly 
     * after \a dependsOn. Use this function to build up the group internal dependencies.
     * \note The node \a dependsOn needs to be already part of the group and not belong to a parent - or
     * subgroup
     *
     * \param callable The callable object that shall be added to this execution group
     * \param dependsOn The node which describes the callable inside the group \node
     * @return void
     */
    Node* createDependendNode( Node* dependsOn, const boost::function<void()>& callable) {

        Node* node = new Node(callable);
        
        setupDependendGroupNode(node, dependsOn);
        return node;
    };
    
        /**
     * @brief Creates a group in this group of connected nodes
     *
     * This function creates a group in the execution group and ensures that it will be executed 
     * after \a this groups setup and before it finishes as a whole. Furthermore it is ensured that
     * \a atStart is executed before and \a atFinish after any other node in the new group.
     *
     * \param atStart The callable object that shall be executed at the new groups startup
     * \param atFinish The callable object that shall be executed when the new group finishs
     * @return The node which describes the group inside the group
     */
    Group* createGroup(const boost::function<void()>& atStart = NULL,
                       const boost::function<void()>& atFinish = NULL) {

        Group* node = new Group(atStart, atFinish);

        setupGroupNode(node);
        return node;
    };

    /**
     * @brief Creates a group in \a this group of connected nodes as dependend on another node
     *
     * This function adds a new group to the execution group and ensures that it will be executed after 
     * \a this groups setup and before it finishes as a whole. It furthermore ensures that it is executed 
     * after the finishing of \a dependsOn. It does however not ensured that it is executed directly 
     * after \a dependsOn. Use this function to build up the group internal dependencies. It is ensured that
     * \a atStart is executed before and \a atFinish after any other node in the new group.
     * 
     * \note The node \a dependsOn needs to be already part of the group and not belong to a parent - or
     * subgroup
     *
     * \param atStart The callable object that shall be executed at the new groups startup
     * \param atFinish The callable object that shall be executed when the new group finishs
     * \param dependsOn The node which describes the callable inside the group \node
     * @return void
     */
    Group* createDependendGroup(Node* dependsOn,
                                const boost::function<void()>& atStart = NULL,
                                const boost::function<void()>& atFinish = NULL ) {

        Group* node = new Group(atStart, atFinish);

        setupDependendGroupNode(node, dependsOn);
        return node;
    };
    
    /**
     * @brief Sets the callable object which is executed at the groups startup
     * 
     * It is ensured that the callable object will be executed before any other group node, however only
     * after all of the groups dependencys are finished.
     * 
     * \param callable The callable object that shall be executed
     * @return void
     */
    void executeAtStart(const boost::function<void()>& callable) {
       m_function = callable;  
    };
    
    /**
     * @brief Sets the callable object which is executed when the group finishs
     * 
     * It is ensured that the callable object will be executed after any other group node, but before any
     * node that dependends on \this will be executed
     * 
     * \param callable The callable object that shall be executed
     * @return void
     */
    void executeAtFinish(const boost::function<void()>& callable) {
        m_final->m_function = callable;
    };

protected:
    Node* m_final;
    
    void setupGroupNode(Node* node) {

        //create a node which is executed after group setup in parallel to all other initial nodes
        //also make sure it is finished before the group is done as a whole by connecting to the
        //final node, which collects all open-end jobs of a group
        m_next.push_back(node);
        node->setupDependendNode(m_final);
    };

    void setupDependendGroupNode(Node* node, Node* dependsOn) {

#ifdef DCM_DEBUG      
        
        //It is only valid to depend on nodes which are inside of this group. If this would not be
        //the case the group execution could be blocked by some unknown dependency which could 
        //make the node execution very unpredictable. Furthermore it violates our guarantees
        dcm_assert(isInGroup(dependsOn));
        
        //The node we depend on must already have a dependend node, eiteher the groups final or at
        //least one other node. If this is not true somehow the setup process has gone wrong.
        dcm_assert(dependsOn->begin() != dependsOn->end());
#endif       
        //we have two possible cases: 
        //1. dependsOn has a connection to m_final which needs to be removed
        //2. dependsOn already has a dependend Node and we can simply add annother one
        if(dependsOn->m_next.front() == m_final) {
            dependsOn->m_next.clear();
            m_final->registeredInput--;
        };
        
        dependsOn->setupDependendNode(node);
        node->setupDependendNode(m_final);
    };

#ifdef DCM_DEBUG
    bool isInGroup(Node* searched) {
        
        BOOST_FOREACH(Node* nn, m_next) {

            if(nn != m_final)
                if(isOrHasAsDependency(nn, searched))
                    return true;
        };
        
        return false;
    };

    bool isOrHasAsDependency(Node* start, Node* searched) {
        
        if(start == searched)
            return true;
        
        bool value = false;        
        BOOST_FOREACH(Node* nn, *start) {
          if(nn == searched)
              return true;
                  
          if(nn != m_final)
            value = isOrHasAsDependency(nn, searched);
          
          if(value)
              return true;
        };
        return false;
    };
#endif
};

//thread pool
class Scheduler {
    
public:
    Scheduler(
#ifdef DCM_USE_MULTITHREADING
        unsigned int n = boost::thread::hardware_concurrency()
#else
        unsigned int n = 0
#endif
    );

    ~Scheduler();
        
    void execute(Node* node);
    void join();

private:
    std::deque< Node* > m_tasks;
#ifdef DCM_USE_MULTITHREADING
    boost::thread_group m_threads;
    boost::mutex m_mutex;
    boost::condition_variable m_cvTask;
    boost::condition_variable m_cvFinished;
    boost::atomic<uint> m_busy;
    bool stop;
#endif

    void enqueue(Node* n);
    void loop();
};

Scheduler::Scheduler(unsigned int n) 
#ifdef DCM_USE_MULTITHREADING
: m_busy(0), stop() 
#endif
{
#ifdef DCM_USE_MULTITHREADING
    for (unsigned int i=0; i<n; ++i)
        m_threads.create_thread(boost::bind(&Scheduler::loop, this));
#endif
}

Scheduler::~Scheduler() {
    
#ifdef DCM_USE_MULTITHREADING
    // set stop-condition
    boost::unique_lock<boost::mutex> lock(m_mutex);
    stop = true;
    m_cvTask.notify_all();
    lock.unlock();

    // all threads terminate, then we're done.
    m_threads.join_all();
#endif
}

void Scheduler::loop() {
    
    while (true) {
        
#ifdef DCM_USE_MULTITHREADING        
        boost::unique_lock<boost::mutex> lock(m_mutex);
        while(!(stop || !m_tasks.empty()))
            m_cvTask.wait(lock);
#endif
        if (!m_tasks.empty())  {

#ifdef DCM_USE_MULTITHREADING 
            // got work. set m_busy.
            m_busy.fetch_add(1, boost::memory_order_relaxed);
#endif
            // pull from queue
            Node* task = m_tasks.front();
            m_tasks.pop_front();
            
#ifdef DCM_USE_MULTITHREADING
            // release lock. run async
            lock.unlock();
#endif           
            task->execute();

            //add the dependend jobs to the queue if they are ready for execution
            BOOST_FOREACH(Node* node, *task) {
                //we need to make the condition >= to ensure nodes without registered inputs are processed
#ifdef DCM_USE_MULTITHREADING
                if(node->handledInputs.fetch_add(1, boost::memory_order_relaxed) >= (node->registeredInput-1)) {
                    node->handledInputs.store(0, boost::memory_order_relaxed);
                    enqueue(node);
                }
#else
                if(++node->handledInputs >= node->registeredInput) {
                    node->handledInputs = 0;
                    enqueue(node);
                }
#endif
            }
#ifdef DCM_USE_MULTITHREADING
            m_busy.fetch_sub(1, boost::memory_order_relaxed);
            m_cvFinished.notify_one();
#endif
        }
#ifdef DCM_USE_MULTITHREADING
        else if (stop)
#else
        else
#endif
            break;
    }

}

void Scheduler::execute(Node* node) {

#ifdef DCM_USE_MULTITHREADING
    boost::unique_lock<boost::mutex> lock(m_mutex);
    m_tasks.clear();
    m_tasks.push_back(node);
    m_cvTask.notify_one();
#else
    m_tasks.push_back(node);
    loop();
#endif
}


// generic function push
void Scheduler::enqueue(Node* node) {
    
#ifdef DCM_USE_MULTITHREADING    
    boost::unique_lock<boost::mutex> lock(m_mutex);
    m_tasks.push_back(node);
    m_cvTask.notify_one();
#else
    m_tasks.push_back(node);
#endif
}

// waits until the queue is empty and all jobs are finished.
void Scheduler::join() {
#ifdef DCM_USE_MULTITHREADING    
    boost::unique_lock<boost::mutex> lock(m_mutex);
    while(!m_tasks.empty() || (m_busy != 0))
        m_cvFinished.wait(lock);
#endif
}

} //shedule
} //details
} //dcm

/**@}*/ //Schedule
/**@}*/ //Core

#endif //DCM_SHEDULER_HPP
