//
//  AdjacencyList.h
//  richCluster
//
//  Created by Sarah on 6/2/25.
//

#ifndef AdjacencyList_h
#define AdjacencyList_h

#include <unordered_map>
#include <unordered_set>

class AdjacencyList {
public:
  AdjacencyList(int n_terms) {
    for (int i=0; i<n_terms; ++i){
      adjList[i] = std::unordered_set<int>();
    }
  }
  
  void addNeighbor(int node, int neighbor);
  bool hasNeighbor(int node, int neighbor);
  
  // return reference to adjList (used to initialize ClusterList)
  const std::unordered_map<int, std::unordered_set<int>>& getAdjList() const {
    return adjList;
  }
  
  size_t size() const {
    return adjList.size();
  }
  
private:
  std::unordered_map<int, std::unordered_set<int>> adjList;
};

#endif /* AdjacencyList_h */
