#ifndef Matcher_H
#define Matcher_H

#include "Graph.h"
#include "Vertex.h"
#include <list>
#include <unordered_map>
#include "Enums.h"
#include "DSU.h"
#include "Blossom.h"
#include "Stack.h"
#include "Frontier.h"
#include "Statistics.h"
#include "StackPusher.h"
#include "ARRStackPusher.h"
#include "concurrentqueue.h"

class Matcher {
public:
    template <typename IT, typename VT, template <typename> class StackType = Stack>
    static void match(Graph<IT, VT>& graph);
    template <typename IT, typename VT, template <typename> class StackType = Stack>
    static void match_parallel(Graph<IT, VT>& graph);
    template <typename IT, typename VT, template <typename> class StackType = Stack>
    static void match_serial(Graph<IT, VT>& graph);
    template <typename IT, typename VT, template <typename> class StackType = Stack>
    static void match(Graph<IT, VT>& graph, Statistics<IT>& stats);
    static void hello_world(int tid);
    template <typename IT, typename VT, template <typename, template <typename> class> class FrontierType, template <typename> class StackType = Stack>
    static void search_master(Graph<IT, VT>& graph, 
                    const size_t V_index,
                    std::vector<FrontierType<IT, StackType>*> & frontiers,
                    bool &foundPath,
                    bool &finished,
                    std::vector<size_t> &read_messages,
                    int tid);
    template <typename IT, typename VT, template <typename, template <typename> class> class FrontierType, template <typename> class StackType = Stack>
    static Vertex<IT> * search(Graph<IT, VT>& graph, 
                    const size_t V_index,
                    FrontierType<IT, StackType> & f);
    template <typename IT, typename VT, template <typename> class StackType>
    static void search_worker(
                    moodycamel::ConcurrentQueue<IT> &worklist,
                    Graph<IT, VT>& graph, 
                    bool &ready,
                    bool &processed,
                    bool &finished,
                    std::vector<size_t> &read_messages,
                    int tid,
                    std::mutex & mtx,
                    std::condition_variable & cv);
    template <typename IT, typename VT, template <typename> class StackType>
    static void search_worker_test(
                    moodycamel::ConcurrentQueue<IT> &worklist,
                    Graph<IT, VT>& graph, 
                    bool &ready,
                    bool &processed,
                    bool &finished,
                    std::vector<size_t> &read_messages,
                    int tid,
                    std::mutex & mtx,
                    std::condition_variable & cv);
    private:
    template <typename IT, typename VT, template <typename> class StackType>
    static void match_master(std::vector<std::thread> &threads,
                                                        unsigned num_threads,
                                                        std::vector<size_t> &read_messages,
                                                        Graph<IT, VT> &graph,
                                                        std::vector<Frontier<IT, StackType>*> &frontiers,
                                                        bool &foundPath,
                                                        bool &finished);

    template <typename IT, typename VT, template <typename, template <typename> class> class FrontierType, template <typename> class StackType = Stack>
    static void augment(Graph<IT, VT>& graph, 
                    Vertex<IT> * TailOfAugmentingPath,
                    FrontierType<IT, StackType> & f);
    template <typename IT, typename VT>
    static void pathThroughBlossom(Graph<IT, VT>& graph, 
                        // V
                        const Vertex<IT> * TailOfAugmentingPath,
                        const Vertex<IT> * TailOfAugmentingPathBase,
                        std::vector<Vertex<IT>> & vertexVector,
                        //std::list<IT> & path,
                        Stack<IT> & path);
};

void Matcher::hello_world(int tid) {
    // Get the Thread ID (TID)
    std::cout << "Hello World from Thread " << tid << std::endl;
}

template <typename IT, typename VT, template <typename> class StackType>
void Matcher::match_master(std::vector<std::thread> &threads,
                                                    unsigned num_threads,
                                                    std::vector<size_t> &read_messages,
                                                    Graph<IT, VT> &graph,
                                                    std::vector<Frontier<IT, StackType>*> &frontiers,
                                                    bool &foundPath,
                                                    bool &finished){
    IT tid = 0;
    auto allocate_start = high_resolution_clock::now();
    //Frontier<IT> f(graph.getN(),graph.getM());
    Frontier<IT,StackType> f(graph.getN(),graph.getM());
    frontiers[tid] = &f;
    auto allocate_end = high_resolution_clock::now();
    auto duration_alloc = duration_cast<milliseconds>(allocate_end - allocate_start);
    std::cout << "Frontier (9|V|+|E|) memory allocation time: "<< duration_alloc.count() << " milliseconds" << '\n';
    Vertex<IT> * TailOfAugmentingPath;
    // Access the graph elements as needed
    for (std::size_t i = 0; i < graph.getN(); ++i) {
        if (graph.matching[i] < 0) {
            //printf("SEARCHING FROM %ld!\n",i);
            // Your matching logic goes here...
            //TailOfAugmentingPath=search(graph,i,f);
            search_master(graph,i,frontiers,foundPath,finished,read_messages,0);
        }
    }
}

template <typename IT, typename VT, template <typename> class StackType>
void Matcher::match(Graph<IT, VT>& graph) {
    auto allocate_start = high_resolution_clock::now();
    //Frontier<IT> f(graph.getN(),graph.getM());
    Frontier<IT,StackType> f(graph.getN(),graph.getM());
    auto allocate_end = high_resolution_clock::now();
    auto duration_alloc = duration_cast<milliseconds>(allocate_end - allocate_start);
    std::cout << "Frontier (9|V|+|E|) memory allocation time: "<< duration_alloc.count() << " milliseconds" << '\n';
    Vertex<IT> * TailOfAugmentingPath;
    // Access the graph elements as needed
    for (std::size_t i = 0; i < graph.getN(); ++i) {
        if (graph.matching[i] < 0) {
            //printf("SEARCHING FROM %ld!\n",i);
            // Your matching logic goes here...
            TailOfAugmentingPath=search(graph,i,f);
            // If not a nullptr, I found an AP.
            if (TailOfAugmentingPath){
                augment(graph,TailOfAugmentingPath,f);
                f.reinit();
                f.clear();
                //printf("FOUND AP!\n");
            } else {
                f.clear();
                //printf("DIDNT FOUND AP!\n");
            }
        }
    }
}


template <typename IT, typename VT, template <typename> class StackType>
void Matcher::match(Graph<IT, VT>& graph, Statistics<IT>& stats) {
    auto allocate_start = high_resolution_clock::now();
    //Frontier<IT> f(graph.getN(),graph.getM());
    Frontier<IT,StackType> f(graph.getN(),graph.getM());
    auto allocate_end = high_resolution_clock::now();
    auto duration_alloc = duration_cast<milliseconds>(allocate_end - allocate_start);
    std::cout << "Frontier (9|V|+|E|) memory allocation time: "<< duration_alloc.count() << " milliseconds" << '\n';
    Vertex<IT> * TailOfAugmentingPath;
    // Access the graph elements as needed
    for (std::size_t i = 0; i < graph.getN(); ++i) {
        if (graph.matching[i] < 0) {
            //printf("SEARCHING FROM %ld!\n",i);
            // Your matching logic goes here...
            auto search_start = high_resolution_clock::now();
            TailOfAugmentingPath=search(graph,i,f);
            auto search_end = high_resolution_clock::now();
            // If not a nullptr, I found an AP.
            if (TailOfAugmentingPath){
                augment(graph,TailOfAugmentingPath,f);
                stats.write_entry(f.path.size() ? (2*f.path.size()-1):0,f.tree.size(),duration_cast<microseconds>(search_end - search_start));
                f.reinit();
                f.clear();
                //printf("FOUND AP!\n");
            } else {
                stats.write_entry(f.path.size() ? (2*f.path.size()-1):0,f.tree.size(),duration_cast<microseconds>(search_end - search_start));
                f.clear();
                //printf("DIDNT FOUND AP!\n");
            }
        }
    }
}


template <typename IT, typename VT, template <typename> class StackType>
void Matcher::search_worker(
                    moodycamel::ConcurrentQueue<IT> &worklist,
                    Graph<IT, VT>& graph, 
                    bool &ready,
                    bool &processed,
                    bool &finished,
                    std::vector<size_t> &read_messages,
                    int tid,
                    std::mutex & mtx,
                    std::condition_variable & cv) {
    Vertex<int64_t> *FromBase,*ToBase, *nextVertex;
    int64_t FromBaseVertexID,ToBaseVertexID;
    IT stackEdge, matchedEdge;
    IT nextVertexIndex;
    IT time;

    //std::cout << "Hello World from Thread " << tid << std::endl;
    auto allocate_start = high_resolution_clock::now();
    //Frontier<IT> f(graph.getN(),graph.getM());
    Frontier<IT,StackType> f(graph.getN(),graph.getM());
    auto allocate_end = high_resolution_clock::now();
    auto duration_alloc = duration_cast<milliseconds>(allocate_end - allocate_start);
    //std::cout << "Thread " << tid << " Frontier (9|V|+|E|) memory allocation time: "<< duration_alloc.count() << " milliseconds" << '\n';


    StackType<IT> &stack = f.stack;
    Stack<IT> &tree = f.tree;
    DisjointSetUnion<IT> &dsu = f.dsu;
    std::vector<Vertex<IT>> & vertexVector = f.vertexVector;
    IT V_index;
    while (!finished) {
        if(!worklist.try_dequeue(V_index))
            continue;

        //std::cout << "Worker thread start" << std::endl;
        std::unique_lock lk(mtx);
        cv.wait(lk, [&] { return ready; });
        //std::cout << "Worker thread " << tid << " is processing data " << V_index << std::endl;
        read_messages[tid]++;
        time = 0;
        nextVertex = &vertexVector[V_index];
        tree.push_back(V_index);
        nextVertex->AgeField=time++;
        StackPusher<IT,VT,StackType>::pushEdgesOntoStack(graph,vertexVector,V_index,stack);
        while(!stack.empty()){
            stackEdge = stack.back();
            stack.pop_back();
            //stack.pop_back();
            // Necessary because vertices dont know their own index.
            // It simplifies vector creation..
            FromBaseVertexID = dsu[Graph<IT,VT>::EdgeFrom(graph,stackEdge)];
            FromBase = &vertexVector[FromBaseVertexID];

            // Necessary because vertices dont know their own index.
            // It simplifies vector creation..
            ToBaseVertexID = dsu[Graph<IT,VT>::EdgeTo(graph,stackEdge)];
            ToBase = &vertexVector[ToBaseVertexID];

            // Edge is between two vertices in the same blossom, continue.
            if (FromBase == ToBase)
                continue;
            if(!FromBase->IsEven()){
                std::swap(FromBase,ToBase);
                std::swap(FromBaseVertexID,ToBaseVertexID);
            }
            // An unreached, unmatched vertex is found, AN AUGMENTING PATH!
            if (!ToBase->IsReached() && !graph.IsMatched(ToBaseVertexID)){
                ToBase->TreeField=stackEdge;
                ToBase->AgeField=time++;
                tree.push_back(ToBaseVertexID);
                //graph.SetMatchField(ToBaseVertexID,stackEdge);
                // I'll let the augment path method recover the path.
                augment(graph,ToBase,f);
                //foundPath = true;
                break;
            } else if (!ToBase->IsReached() && graph.IsMatched(ToBaseVertexID)){
                ToBase->TreeField=stackEdge;
                ToBase->AgeField=time++;
                tree.push_back(ToBaseVertexID);

                matchedEdge=graph.GetMatchField(ToBaseVertexID);
                nextVertexIndex = Graph<IT,VT>::Other(graph,matchedEdge,ToBaseVertexID);
                nextVertex = &vertexVector[nextVertexIndex];
                nextVertex->AgeField=time++;
                tree.push_back(nextVertexIndex);

                //Graph<IT,VT>::pushEdgesOntoStack(graph,vertexVector,nextVertexIndex,stack,matchedEdge);
                StackPusher<IT,VT,StackType>::pushEdgesOntoStack(graph,vertexVector,nextVertexIndex,stack,matchedEdge);

            } else if (ToBase->IsEven()) {
                // Shrink Blossoms
                // Not sure if this is wrong or the augment method is wrong
                Blossom::Shrink(graph,stackEdge,dsu,vertexVector,stack);
            }
        }
        f.reinit();
        f.clear();

        // Send data back to master thread
        processed = true;
        //std::cout << "Worker thread signals data processing completed" << std::endl;

        // Manual unlocking is done before notifying, to avoid waking up
        // the waiting thread only to block again (see notify_one for details)
        lk.unlock();
        // The worker thread has done the work,
        // Notify the master thread to continue the work.
        cv.notify_one();
    }
}


template <typename IT, typename VT, template <typename> class StackType>
void Matcher::search_worker_test(
                    moodycamel::ConcurrentQueue<IT> &worklist,
                    Graph<IT, VT>& graph, 
                    bool &ready,
                    bool &processed,
                    bool &finished,
                    std::vector<size_t> &read_messages,
                    int tid,
                    std::mutex & mtx,
                    std::condition_variable & cv) {
    Vertex<int64_t> *FromBase,*ToBase, *nextVertex;
    int64_t FromBaseVertexID,ToBaseVertexID;
    IT stackEdge, matchedEdge;
    IT nextVertexIndex;
    IT time;

    //std::cout << "Hello World from Thread " << tid << std::endl;
    auto allocate_start = high_resolution_clock::now();
    //Frontier<IT> f(graph.getN(),graph.getM());
    Frontier<IT,StackType> f(graph.getN(),graph.getM());
    auto allocate_end = high_resolution_clock::now();
    auto duration_alloc = duration_cast<milliseconds>(allocate_end - allocate_start);
    //std::cout << "Thread " << tid << " Frontier (9|V|+|E|) memory allocation time: "<< duration_alloc.count() << " milliseconds" << '\n';


    StackType<IT> &stack = f.stack;
    Stack<IT> &tree = f.tree;
    DisjointSetUnion<IT> &dsu = f.dsu;
    std::vector<Vertex<IT>> & vertexVector = f.vertexVector;
    IT V_index;
    Vertex<IT> * TailOfAugmentingPath;
    // Access the graph elements as needed
    for (std::size_t i = 0; i < graph.getN(); ++i) {
        if (graph.matching[i] < 0) {
            //printf("SEARCHING FROM %ld!\n",i);
            // Your matching logic goes here...
            TailOfAugmentingPath=search(graph,i,f);
            // If not a nullptr, I found an AP.
            if (TailOfAugmentingPath){
                augment(graph,TailOfAugmentingPath,f);
                f.reinit();
                f.clear();
                //printf("FOUND AP!\n");
            } else {
                f.clear();
                //printf("DIDNT FOUND AP!\n");
            }
        }
    }
    finished = true;
}


template <typename IT, typename VT, template <typename, template <typename> class> class FrontierType, template <typename> class StackType = Stack>
void Matcher::search_master(Graph<IT, VT>& graph, 
                    const size_t V_index,
                    std::vector<FrontierType<IT, StackType>*> & frontiers,
                    bool &foundPath,
                    bool &finished,
                    std::vector<size_t> &read_messages,
                    int tid){
    Vertex<int64_t> *FromBase,*ToBase, *nextVertex;
    int64_t FromBaseVertexID,ToBaseVertexID;
    IT stackEdge, matchedEdge;
    IT nextVertexIndex;
    IT time = 0;
    FrontierType<IT, StackType> & frontier = *(frontiers[tid]);
    StackType<IT> &stack = frontier.stack;
    Stack<IT> &tree = frontier.tree;
    DisjointSetUnion<IT> &dsu = frontier.dsu;
    std::vector<Vertex<IT>> & vertexVector = frontier.vertexVector;
    nextVertex = &vertexVector[V_index];
    tree.push_back(V_index);
    nextVertex->AgeField=time++;
    ARRStackPusher::pushEdgesOntoStack<IT,VT,FrontierType>(graph,vertexVector,V_index,frontiers,tid);
    StackPusher<IT,VT,StackType>::pushEdgesOntoStack(graph,vertexVector,V_index,stack);
    while(!stack.empty()){
        stackEdge = stack.back();
        stack.pop_back();
        //stack.pop_back();
        // Necessary because vertices dont know their own index.
        // It simplifies vector creation..
        FromBaseVertexID = dsu[Graph<IT,VT>::EdgeFrom(graph,stackEdge)];
        FromBase = &vertexVector[FromBaseVertexID];

        // Necessary because vertices dont know their own index.
        // It simplifies vector creation..
        ToBaseVertexID = dsu[Graph<IT,VT>::EdgeTo(graph,stackEdge)];
        ToBase = &vertexVector[ToBaseVertexID];

        // Edge is between two vertices in the same blossom, continue.
        if (FromBase == ToBase)
            continue;
        if(!FromBase->IsEven()){
            std::swap(FromBase,ToBase);
            std::swap(FromBaseVertexID,ToBaseVertexID);
        }
        // An unreached, unmatched vertex is found, AN AUGMENTING PATH!
        if (!ToBase->IsReached() && !graph.IsMatched(ToBaseVertexID)){
            ToBase->TreeField=stackEdge;
            ToBase->AgeField=time++;
            tree.push_back(ToBaseVertexID);
            //graph.SetMatchField(ToBaseVertexID,stackEdge);
            // I'll let the augment path method recover the path.
            augment(graph,ToBase,frontier);
            frontier.reinit();
            frontier.clear();
            return;
        } else if (!ToBase->IsReached() && graph.IsMatched(ToBaseVertexID)){
            ToBase->TreeField=stackEdge;
            ToBase->AgeField=time++;
            tree.push_back(ToBaseVertexID);

            matchedEdge=graph.GetMatchField(ToBaseVertexID);
            nextVertexIndex = Graph<IT,VT>::Other(graph,matchedEdge,ToBaseVertexID);
            nextVertex = &vertexVector[nextVertexIndex];
            nextVertex->AgeField=time++;
            tree.push_back(nextVertexIndex);

            //Graph<IT,VT>::pushEdgesOntoStack(graph,vertexVector,nextVertexIndex,stack,matchedEdge);
            StackPusher<IT,VT,StackType>::pushEdgesOntoStack(graph,vertexVector,nextVertexIndex,stack,matchedEdge);

        } else if (ToBase->IsEven()) {
            // Shrink Blossoms
            // Not sure if this is wrong or the augment method is wrong
            Blossom::Shrink(graph,stackEdge,dsu,vertexVector,stack);
        }
    }
    frontier.clear();
    return;
}

template <typename IT, typename VT, template <typename, template <typename> class> class FrontierType, template <typename> class StackType = Stack>
Vertex<IT> * Matcher::search(Graph<IT, VT>& graph, 
                    const size_t V_index,
                    FrontierType<IT, StackType> & f) {
    Vertex<int64_t> *FromBase,*ToBase, *nextVertex;
    int64_t FromBaseVertexID,ToBaseVertexID;
    IT stackEdge, matchedEdge;
    IT nextVertexIndex;
    IT time = 0;
    StackType<IT> &stack = f.stack;
    Stack<IT> &tree = f.tree;
    DisjointSetUnion<IT> &dsu = f.dsu;
    std::vector<Vertex<IT>> & vertexVector = f.vertexVector;
    //auto inserted = vertexMap.try_emplace(V_index,Vertex<IT>(time++,Label::EvenLabel));
    nextVertex = &vertexVector[V_index];
    tree.push_back(V_index);
    nextVertex->AgeField=time++;
    // Push edges onto stack, breaking if that stackEdge is a solution.
    //Graph<IT,VT>::pushEdgesOntoStack(graph,vertexVector,V_index,stack);
    StackPusher<IT,VT,StackType>::pushEdgesOntoStack(graph,vertexVector,V_index,stack);
    while(!stack.empty()){
        stackEdge = stack.back();
        stack.pop_back();
        //stack.pop_back();
        // Necessary because vertices dont know their own index.
        // It simplifies vector creation..
        FromBaseVertexID = dsu[Graph<IT,VT>::EdgeFrom(graph,stackEdge)];
        FromBase = &vertexVector[FromBaseVertexID];

        // Necessary because vertices dont know their own index.
        // It simplifies vector creation..
        ToBaseVertexID = dsu[Graph<IT,VT>::EdgeTo(graph,stackEdge)];
        ToBase = &vertexVector[ToBaseVertexID];

        // Edge is between two vertices in the same blossom, continue.
        if (FromBase == ToBase)
            continue;
        if(!FromBase->IsEven()){
            std::swap(FromBase,ToBase);
            std::swap(FromBaseVertexID,ToBaseVertexID);
        }
        // An unreached, unmatched vertex is found, AN AUGMENTING PATH!
        if (!ToBase->IsReached() && !graph.IsMatched(ToBaseVertexID)){
            ToBase->TreeField=stackEdge;
            ToBase->AgeField=time++;
            tree.push_back(ToBaseVertexID);
            //graph.SetMatchField(ToBaseVertexID,stackEdge);
            // I'll let the augment path method recover the path.
            return ToBase;
        } else if (!ToBase->IsReached() && graph.IsMatched(ToBaseVertexID)){
            ToBase->TreeField=stackEdge;
            ToBase->AgeField=time++;
            tree.push_back(ToBaseVertexID);

            matchedEdge=graph.GetMatchField(ToBaseVertexID);
            nextVertexIndex = Graph<IT,VT>::Other(graph,matchedEdge,ToBaseVertexID);
            nextVertex = &vertexVector[nextVertexIndex];
            nextVertex->AgeField=time++;
            tree.push_back(nextVertexIndex);

            //Graph<IT,VT>::pushEdgesOntoStack(graph,vertexVector,nextVertexIndex,stack,matchedEdge);
            StackPusher<IT,VT,StackType>::pushEdgesOntoStack(graph,vertexVector,nextVertexIndex,stack,matchedEdge);

        } else if (ToBase->IsEven()) {
            // Shrink Blossoms
            // Not sure if this is wrong or the augment method is wrong
            Blossom::Shrink(graph,stackEdge,dsu,vertexVector,stack);
        }
    }
    return nullptr;
}

template <typename IT, typename VT, template <typename, template <typename> class> class FrontierType, template <typename> class StackType = Stack>
void Matcher::augment(Graph<IT, VT>& graph, 
                    Vertex<IT> * TailOfAugmentingPath,
                    FrontierType<IT, StackType> & f) {

    DisjointSetUnion<IT> &dsu = f.dsu;
    std::vector<Vertex<IT>> & vertexVector = f.vertexVector;
    //std::list<IT> path;
    Stack<IT> & path = f.path;
    IT edge;
    // W
    Vertex<IT>*nextVertex;
    Vertex<IT>*nextVertexBase;
    do
    {
        //ListPut(Tree(V), P);
        edge = TailOfAugmentingPath->TreeField;
        path.push_back(edge);

        //W = Other(Tree(V), V);
        ptrdiff_t TailOfAugmentingPath_VertexID = TailOfAugmentingPath - &vertexVector[0];
        auto nextVertexID = Graph<IT,VT>::Other(graph,edge,TailOfAugmentingPath_VertexID);
        nextVertex = &vertexVector[nextVertexID];

        //B = Base(Blossom(W));
        auto nextVertexBaseID = dsu[nextVertexID];  
        nextVertexBase = &vertexVector[nextVertexBaseID];
        
        // Path(W, B, P);
        pathThroughBlossom(graph,nextVertex,nextVertexBase,vertexVector,path);

        //V = Other(Match(B), B);
        ptrdiff_t nextVertexBase_VertexID = nextVertexBase - &vertexVector[0];
        if (graph.IsMatched(nextVertexBase_VertexID))
            TailOfAugmentingPath = &vertexVector[Graph<IT,VT>::Other(graph,graph.GetMatchField(nextVertexBase_VertexID),nextVertexBase_VertexID)];
        else 
            TailOfAugmentingPath = nullptr;
    } while (TailOfAugmentingPath != nullptr);
    // Print the list of integers
    for (auto E : path) {
        //Match(EdgeFrom(E)) = E;
        graph.SetMatchField(Graph<IT,VT>::EdgeFrom(graph,E),E);
        //Match(EdgeTo(E)) = E;
        graph.SetMatchField(Graph<IT,VT>::EdgeTo(graph,E),E);
    }
}


template <typename IT, typename VT>
void Matcher::pathThroughBlossom(Graph<IT, VT>& graph, 
                    // V
                    const Vertex<IT> * TailOfAugmentingPath,
                    const Vertex<IT> * TailOfAugmentingPathBase,
                    std::vector<Vertex<IT>> & vertexVector,
                    //std::list<IT> & path,
                    Stack<IT> & path) {
    // W
    Vertex<IT>*nextVertex;
    // if (V != B)
    if (TailOfAugmentingPath != TailOfAugmentingPathBase)
    {
        if (TailOfAugmentingPath->IsOdd())
        {
            // Path(Shore(V), Other(Match(V), V), P);
            ptrdiff_t TailOfAugmentingPath_VertexID = TailOfAugmentingPath - &vertexVector[0];
            pathThroughBlossom(graph,
                                &vertexVector[TailOfAugmentingPath->ShoreField],
                                &vertexVector[Graph<IT,VT>::Other(graph,graph.GetMatchField(TailOfAugmentingPath_VertexID),TailOfAugmentingPath_VertexID)],
                                vertexVector,
                                path);
            //ListPut(Bridge(V), P);
            path.push_back(TailOfAugmentingPath->BridgeField);
            
            //Path(Other(Bridge(V), Shore(V)), B, P);
            pathThroughBlossom(graph,
                                &vertexVector[Graph<IT,VT>::Other(graph,TailOfAugmentingPath->BridgeField,TailOfAugmentingPath->ShoreField)],
                                TailOfAugmentingPathBase,
                                vertexVector,
                                path);
        }
        else if (TailOfAugmentingPath->IsEven())
        {
            //W = Other(Match(V), V);
            ptrdiff_t TailOfAugmentingPath_VertexID = TailOfAugmentingPath - &vertexVector[0];
            nextVertex=&vertexVector[Graph<IT,VT>::Other(graph,graph.GetMatchField(TailOfAugmentingPath_VertexID),TailOfAugmentingPath_VertexID)];
            
            //ListPut(Tree(W), P);
            path.push_back(nextVertex->TreeField);

            //Path(Other(Tree(W), W), B, P);
            ptrdiff_t nextVertex_VertexID = nextVertex - &vertexVector[0];
            pathThroughBlossom(graph,
                                &vertexVector[Graph<IT,VT>::Other(graph,nextVertex->TreeField,nextVertex_VertexID)],
                                TailOfAugmentingPathBase,
                                vertexVector,
                                path);
        }
        else{
            ptrdiff_t TailOfAugmentingPath_VertexID = TailOfAugmentingPath - &vertexVector[0];
            ptrdiff_t TailOfAugmentingPathBase_VertexID = TailOfAugmentingPathBase - &vertexVector[0];
            std::cerr << "(Path) Internal error. TailOfAugmentingPath_VertexID: " << TailOfAugmentingPath_VertexID<< " TailOfAugmentingPathBase_VertexID: " << TailOfAugmentingPathBase_VertexID << std::endl;
            exit(1);
        }
    }
}
// Needs to be after Matcher definition.
// Not sure how to fix this.
#include "ThreadFactory2.h"
// Mutex/CV inspired from https://leimao.github.io/blog/CPP-Condition-Variable/
template <typename IT, typename VT, template <typename> class StackType>
void Matcher::match_parallel(Graph<IT, VT>& graph) {
    bool finished = false;
    
    bool ready = false;
    bool processed = false;
    bool foundPath = false;

    // This will take the place of foundPath/finished on a vertex
    std::mutex mtx;
    std::condition_variable cv;
    // Multi-producer; Multi-consumer queue, aka worklist
    // Start off by using it as a dynamic allocator of work.
    size_t capacity = 1;
    moodycamel::ConcurrentQueue<IT> worklist{capacity};
    // 8 workers.
    constexpr unsigned num_threads = 1;
    std::vector<std::thread> workers(num_threads);
    std::vector<size_t> read_messages;
    read_messages.resize(num_threads);
    // Access the graph elements as needed
    ThreadFactory::create_threads_concurrentqueue_baseline<IT,VT,StackType>(workers, num_threads,read_messages,worklist,graph,ready,processed,finished,mtx,cv);

    IT written_messages = 0;
    cpu_set_t my_set;
    CPU_ZERO(&my_set);
    CPU_SET(0, &my_set);
    if (sched_setaffinity(0, sizeof(cpu_set_t), &my_set)) {
        std::cout << "sched_setaffinity error: " << strerror(errno) << std::endl;
    }
    auto match_start = high_resolution_clock::now();


    for (std::size_t i = 0; i < graph.getN(); ++i) {
        if (graph.matching[i] < 0) {
            ready = false;
            processed = false;
            //std::cout << "Master thread start " << i << std::endl;
            // Send data to the worker thread.
            {
                worklist.enqueue(i);
                std::lock_guard lk(mtx);
                ready = true;
                //std::cout << "Master thread signals data ready for processing"
                //        << std::endl;
            }
            // The mater thread has done the preliminary work,
            // Notify the worker thread to continue the work.
            cv.notify_one();

            // Wait for the worker.
            {
                std::unique_lock lk(mtx);
                cv.wait(lk, [&] { return processed; });
            }
            //std::cout << "Back in master thread, data = " << i << std::endl;
        }
    }

    //Matcher::match_master<IT, VT, StackType>(workers, num_threads,read_messages,graph,frontiers,foundPath,finished);
    auto match_end = high_resolution_clock::now();

    auto duration = duration_cast<milliseconds>(match_end - match_start);

    finished=true;

    
    print_results(BenchResult{num_threads, written_messages, read_messages, duration});

    for (auto& t : workers) {
        t.join();
    }
    
    return;
}

// Mutex/CV inspired from https://leimao.github.io/blog/CPP-Condition-Variable/
template <typename IT, typename VT, template <typename> class StackType>
void Matcher::match_serial(Graph<IT, VT>& graph) {
    bool finished = false;
    
    bool ready = false;
    bool processed = false;
    bool foundPath = false;

    // This will take the place of foundPath/finished on a vertex
    std::mutex mtx;
    std::condition_variable cv;
    // Multi-producer; Multi-consumer queue, aka worklist
    // Start off by using it as a dynamic allocator of work.
    size_t capacity = 1;
    moodycamel::ConcurrentQueue<IT> worklist{capacity};
    // 8 workers.
    constexpr unsigned num_threads = 1;
    std::vector<std::thread> workers(num_threads);
    std::vector<size_t> read_messages;
    read_messages.resize(num_threads);
    // Access the graph elements as needed
    ThreadFactory::create_threads_concurrentqueue_baseline<IT,VT,StackType>(workers, num_threads,read_messages,worklist,graph,ready,processed,finished,mtx,cv);

    IT written_messages = 0;
    cpu_set_t my_set;
    CPU_ZERO(&my_set);
    CPU_SET(0, &my_set);
    if (sched_setaffinity(0, sizeof(cpu_set_t), &my_set)) {
        std::cout << "sched_setaffinity error: " << strerror(errno) << std::endl;
    }
    auto match_start = high_resolution_clock::now();
    while(!finished){}

    //Matcher::match_master<IT, VT, StackType>(workers, num_threads,read_messages,graph,frontiers,foundPath,finished);
    auto match_end = high_resolution_clock::now();

    auto duration = duration_cast<milliseconds>(match_end - match_start);

    finished=true;

    
    print_results(BenchResult{num_threads, written_messages, read_messages, duration});

    for (auto& t : workers) {
        t.join();
    }
    
    return;
}
#endif
