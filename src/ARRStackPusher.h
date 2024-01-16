#ifndef ARR_STACK_PUSHER_H
#define ARR_STACK_PUSHER_H
#include "Graph.h"
#include "Vertex.h"
#include "Stack.h"
#include <vector>
class ARRStackPusher {
public:
    template <typename IT, typename VT, template <typename, template <typename> class> class FrontierType, template <typename> class StackType = Stack>
    static bool pushEdgesOntoStack(const Graph<IT, VT>& graph,
                                    std::vector<Vertex<IT>>& vertexVector,
                                    IT V_index,
                                    std::vector<FrontierType<IT, StackType>*> & frontiers,
                                    int tid,
                                    IT optionalEdge1=-1,
                                    IT optionalEdge2=-1);
};

template <typename IT, typename VT, template <typename, template <typename> class> class FrontierType, template <typename> class StackType = Stack>
bool ARRStackPusher::pushEdgesOntoStack(const Graph<IT, VT>& graph,
                                    std::vector<Vertex<IT>>& vertexVector,
                                    IT V_index,
                                    std::vector<FrontierType<IT, StackType>*> & frontiers,
                                    int tid,
                                    IT optionalEdge1,
                                    IT optionalEdge2){
    IT nextVertexIndex;
    Vertex<IT>* nextVertex;
    // Push edges onto stack, breaking if that stackEdge is a solution.
    for (IT start = graph.indptr[V_index]; start < graph.indptr[V_index + 1]; ++start) {
        // For blossom contraction, need to skip repushing the matched & tree edges
        if (graph.indices[start] == optionalEdge1 || graph.indices[start] == optionalEdge2)
            continue;
        //stack.push_back(graph.indices[start]);

        nextVertexIndex = Graph<IT, VT>::Other(graph, graph.indices[start], V_index);

        nextVertex = &vertexVector[nextVertexIndex];
        if (!nextVertex->IsReached() && !graph.IsMatched(nextVertexIndex))
            return true;
    }
    return false;
}
#endif