//
//  DistanceMetric.h
//  richCluster
//
//  Created by Sarah on 6/2/25.
//

#ifndef DistanceMetric_h
#define DistanceMetric_h

#include <unordered_set>
#include <string>

class DistanceMetric {
public:
  DistanceMetric(std::string distanceMetric, double distanceCutoff)
    : metric(distanceMetric), cutoff(distanceCutoff) {};
  
  double computeDistance(
      const std::unordered_set<std::string>& t1_genes,
      const std::unordered_set<std::string>& t2_genes,
      int totalGeneCount);
  double getCutoff() { return cutoff; };
  
private:
  std::string metric;
  double cutoff;
  
  // methods
  double getKappa(const std::unordered_set<std::string>& t1_genes,
                  const std::unordered_set<std::string>& t2_genes,
                  int totalGeneCount);
  double getJaccard(const std::unordered_set<std::string>& t1_genes,
                    const std::unordered_set<std::string>& t2_genes);
};

#endif /* DistanceMetric_h */
