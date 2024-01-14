#ifndef FRONTIER_H
#define FRONTIER_H
#include <vector>
#include "Vertex.h"
#include "Stack.h"
#include "DSU.h"

template <typename IT>
class Frontier  {
public:
    Frontier(size_t N, size_t M);
    void reinit();
    void clear();

    // Other member functions...
    std::vector<Vertex<IT>> vertexVector;
    Stack<IT> stack;
    Stack<IT> tree;
    Stack<IT> path;
    DisjointSetUnion<IT> dsu;
};

// Constructor
template <typename IT>
Frontier<IT>::Frontier(size_t N, size_t M): vertexVector(N), tree(N), path(M), stack(M){
    dsu.reset(N);
}

// Constructor
template <typename IT>
void Frontier<IT>::reinit(){       
    for (auto V : tree) {
        vertexVector[V].TreeField=-1;
        vertexVector[V].BridgeField=-1;
        vertexVector[V].ShoreField=-1;
        vertexVector[V].AgeField=-1;
        dsu.link[V]=V;
        dsu.directParent[V]=-1;
        dsu.groupRoot[V]=V;
        dsu.size[V]=1;
    }
}


// Constructor
template <typename IT>
void Frontier<IT>::clear(){       
    stack.clear();
    tree.clear();
    path.clear();
}

#endif // FRONTIER_H
