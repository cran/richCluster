//
//  ClusterList.h
//  richCluster
//
//  Created by Sarah on 6/2/25.
//

#ifndef ClusterList_h
#define ClusterList_h

#include <list>
#include <set>
#include <string>
#include <vector>

class ClusterList{
public:
  // C9 (converge ledger, 2026-09-04): std::set, NOT std::unordered_set.  A
  // std::set iterates in ascending order on every standard library, so seed
  // growth, linkage summation, merge choice and export no longer depend on the
  // hash-set traversal that differs between libc++ (macOS) and libstdc++
  // (Linux, mingw).  LinkageMethod::Cluster must stay the same type.
  using Cluster = std::set<int>;
  using ClusterIt = std::list<Cluster>::iterator;
  
  ClusterList(std::vector<std::string>& terms): terms(terms) {};
  void addCluster(Cluster cluster) {clusterList.push_back(cluster);};
  void removeCluster(ClusterIt it) {clusterList.erase(it);};
  void mergeClusters(ClusterIt it1, ClusterIt it2) {
    if (it1 == it2) return;
    it1->insert(it2->begin(), it2->end());
    clusterList.erase(it2);
  }
  std::list<Cluster>& getList() {return clusterList;};
  Rcpp::DataFrame export_r() const;
  void deduplicate();
  size_t size() const { return clusterList.size(); }
  
private:
  std::vector<std::string> terms;
  std::list<Cluster> clusterList;
};

#endif /* ClusterList_h */
