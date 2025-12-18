//
//  LinkageMethod.h
//  richCluster
//
//  Created by Sarah on 6/2/25.
//

#ifndef LinkageMethod_h
#define LinkageMethod_h

#include <unordered_set>
#include <string>
#include <functional>

class LinkageMethod {
public:
  using Cluster = std::unordered_set<int>;
  LinkageMethod(std::string linkageMethod, double linkageCutoff, std::function<double(int, int)> distFn):
    method(linkageMethod), cutoff(linkageCutoff), distFct(distFn) {};
  double computeLinkage(
      const Cluster& cluster1,
      const Cluster& cluster2);
  double getCutoff() { return cutoff; };
  
private:
  std::string method;
  double cutoff;
  std::function<double(int, int)> distFct; // useful reference 
  
  // supported methods
  double single(const Cluster& cluster1,
                const Cluster& cluster2);
  double complete(const Cluster& cluster1,
                  const Cluster& cluster2);
  double average(const Cluster& cluster1,
                 const Cluster& cluster2);
  double ward(const Cluster& cluster1,
              const Cluster& cluster2);
};

#endif /* LinkageMethod_h */
