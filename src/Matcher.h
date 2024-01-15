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

class Matcher {
public:
    template <typename IT, typename VT, template <typename> class StackType = Stack>
    static void match(Graph<IT, VT>& graph);
    template <typename IT, typename VT, template <typename> class StackType = Stack>
    static void match(Graph<IT, VT>& graph, Statistics<IT>& stats);

private:
    template <typename IT, typename VT, template <typename, template <typename> class> class FrontierType, template <typename> class StackType = Stack>
    static Vertex<IT> * search(Graph<IT, VT>& graph, 
                    const size_t V_index,
                    FrontierType<IT, StackType> & f);
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
        /*
        std::optional<IT> stackEdgeOpt = stack.pop();
        if (stackEdgeOpt.has_value())
            stackEdge = stackEdgeOpt.value();
        else
            continue;
        */
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


#endif
