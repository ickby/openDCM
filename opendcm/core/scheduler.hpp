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

#ifdef DCM_USE_MULTITHREADING
#include <tbb/parallel_for_each.h>

#endif

#include "defines.hpp"

#include <tbb/parallel_for.h>
#include <tbb/parallel_for_each.h>
#include <tbb/flow_graph.h>
#include <tbb/tbb.h>


namespace dcm {
namespace shedule {

struct Executable {
    
    virtual ~Executable() {};
    
    void operator()() {execute();};
    virtual void execute() = 0;
};

/**
 * @brief Execute arbitrary functors inside of the sheduler environment
 * 
 * To allow the excecution of arbitrary functors inside the sheduled and parrallel framework this 
 * abstraction class is created. It tores the function object and executes it from within the virtual
 * execute() call. Note that the functor gets copied, hence it must be copy-constructable. Also, if 
 * a state is needed, it must e stored by pointer or reference inside the functor
 */
template<typename T>
struct Functor : public Executable {
    
    Functor(const T& t) : m_functor(t) {};    

    void operator()() {m_functor();};  
    virtual void execute() { operator()();};
    
protected:
    T m_functor;
};

/**
 * @brief Make a Executable from arbitrary functors
 * 
 * To allow the excecution of arbitrary functors this 
 * function can be used. It takes any parameter-free lambda [](){} and returnes an executable, 
 * which can be used within the shedule parallel environment
 */
template<typename T>
std::shared_ptr<Executable> make_executable(const T& functor) {
    return std::make_shared<Functor<T>>(functor);
};

template<typename Vector>
struct _SequentialVector : public Executable {
  
    typedef typename Vector::value_type value_type;
    
    void operator()() {
        for(value_type ex : m_executables) 
            ex->execute();
    };
    
    virtual void execute() {
        operator()();
    };
    
    template<typename T> 
    void add(const T& t) {
        m_executables.push_back(make_executable(t));
    };
    
    void addExecutable(value_type ex) {
        m_executables.push_back(ex);
    };
    
    void append(Vector& vec) {
        m_executables.insert(m_executables.end(), vec.begin(), vec.end());  
    };
    
    int size() { return m_executables.size();};
    bool empty() {return m_executables.empty();};
    
protected:
    Vector m_executables;    
};

/**
 * @brief Sequential execution of multiple Executables
 * 
 * This class stores different executables and processes them sequentially. Note that passed 
 * executable pointers are afterwards owned by the Vector object which delets it when destroyed.
 */
struct SequentialVector : public _SequentialVector<std::vector<std::shared_ptr<Executable>>> {};

/**
 * @brief Sequential execution of Executables with concurrent adding
 * 
 * This class stores different executables and processes them sequentially. It allowes concurrent 
 * adding of executables, which means it can be filled up from the shedule parallel functions.
 * Note that passed executable pointers are afterwards owned by the Vector object which delets it 
 * when destroyed.
 */
struct SequentialConcurrentVector : public _SequentialVector<tbb::concurrent_vector<std::shared_ptr<Executable>>> {};

template<typename Vector>
struct _ParallelVector : public _SequentialVector<Vector> {
    
    void operator()() {
         tbb::parallel_for_each(_SequentialVector<Vector>::m_executables.begin(), _SequentialVector<Vector>::m_executables.end(), 
                            [&](std::shared_ptr<Executable> exe) {
             exe->execute();
         });
    };
    
    virtual void execute() {
        operator()();
    };
};

/**
 * @brief Parallel unordered execution 
 * 
 * This class stores different executables and processes them in parallel. Note that passed 
 * executable pointers are afterwards owned by the Vector object which delets it when destroyed.
 */
struct ParallelVector : public _ParallelVector<std::vector<std::shared_ptr<Executable>>> {};

/**
 * @brief Parallel unordered execution with concurrent adding
 * 
 * This class stores different executables and processes them in parallel. It also allows concurrent
 * adding of executables, which means it can be filled from shedule's parallel functions. 
 * Note that passed executable pointers are afterwards owned by the Vector object which delets it when 
 * destroyed.
 */
struct ParallelConcurrentVector : public _ParallelVector<tbb::concurrent_vector<std::shared_ptr<Executable>>> {};

template<typename Vector>
struct _HugeParallelVector : public _SequentialVector<Vector> {
    
    void operator()() {
         tbb::parallel_for( tbb::blocked_range<int>( 1, _SequentialVector<Vector>::m_executables.size()), 
                            [&](const tbb::blocked_range<int>& range) {
             for(int i=range.begin(); i!=range.end(); ++i)
                 _SequentialVector<Vector>::m_executables[i]->execute();
         });
    };
    
    virtual void execute() {
        operator()();
    };
};

/**
 * @brief Parallel unordered execution of large numbers of execuatbles
 * 
 * This class stores different executables and processes them in parallel. It provides a perforamnce
 * benefit to the \ref ParallelVector in case of many executables (>100).  Note that passed 
 * executable pointers are afterwards owned by the Vector object which delets it when destroyed.
 */
struct HugeParallelVector : public _HugeParallelVector<std::vector<std::shared_ptr<Executable>>> {};

/**
 * @brief Parallel unordered execution of large numbers of execuatbles
 * 
 * This class stores different executables and processes them in parallel. It provides a perforamnce
 * benefit to the \ref ParallelVector in case of many executables (>100). This vector also provides 
 * concurrent adding of executables, it can be filled with the parallel shedule functions.
 * Note that passed executable pointers are afterwards owned by the Vector object which delets it 
 * when destroyed.
 */
struct HugeParallelConcurrentVector : public _HugeParallelVector<tbb::concurrent_vector<std::shared_ptr<Executable>>> {};

//encapsulates a tbb flow graph and is responsible for managing the nodes lifetime
struct FlowGraph : public Executable {
       
    typedef tbb::flow::continue_msg                             ContinueMessage;
    typedef tbb::flow::continue_node< tbb::flow::continue_msg > Node;
    typedef tbb::flow::broadcast_node<tbb::flow::continue_msg>  StartNode;
       
    FlowGraph() : m_graph(new tbb::flow::graph()) {};
    virtual ~FlowGraph() {};
    
    void operator()() {
        m_start.try_put(tbb::flow::continue_msg());
        m_graph->wait_for_all();
    }
    
    virtual void execute() {
        operator()();
    }
        
    template<typename Action>
    Node& newActionNode(Action a) {
        
        m_nodes.emplace_back(Node(*m_graph, a));
        return m_nodes.back();
    };
    
    template<typename Action>
    Node& newInitialActionNode(Action a) {
        
        m_nodes.emplace_back(Node(*m_graph, a));
        connect(m_start, m_nodes.back());
        return m_nodes.back();
    };
    
    StartNode& getBroadcastNode() {        
        return m_start;
    };
    
    template<typename N1, typename N2>
    void connect(N1& n1, N2& n2) {
        tbb::flow::make_edge(n1, n2);
    };

private:
    std::vector<Node>                 m_nodes;
    std::unique_ptr<tbb::flow::graph> m_graph;
    StartNode                         m_start = StartNode(*m_graph);
};
    
//functions for parallel execution

/**
 * @brief Parallel iteration between start and end iterator
 * 
 * This function applies each iterator between start and end to the functor in parallel. This is 
 * intended for use with iteratable containers like std::vector and their begin() and end() functions.
 */
template<typename Iterator, typename Functor>
void for_each(const Iterator& start, const Iterator& end, const Functor& func) {
    tbb::parallel_for_each(start, end, func);
};

/**
 * @brief Parallel iteration over a numerical range
 * 
 * This function applies the functor to each number between start and end. The increment is defined
 * by step. This is equaivalent to the default c++ for loop for(i=start;i<end; i+=step) f(i)
 */
template<typename Type, typename Functor>
void for_each(Type start, Type end, Type step, const Functor& func) {
    tbb::parallel_for(start, end, step, func);
};

} //details
} //dcm

/**@}*/ //Schedule
/**@}*/ //Core

#endif //DCM_SHEDULER_HPP
