//
//  richCluster.h
//  richCluster
//
//  Created by Sarah on 6/1/25.
//

#ifndef richCluster_h
#define richCluster_h

#include <Rcpp.h>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <string>
#include <functional>

#include "DistanceMatrix.h"
#include "AdjacencyList.h"
#include "ClusterList.h"
#include "DistanceMetric.h"
#include "LinkageMethod.h"


class richCluster {
public:
  richCluster(Rcpp::CharacterVector r_terms,
              Rcpp::CharacterVector r_geneIDs,
              std::string distanceMetric, double distanceCutoff,
              std::string linkageMethod, double linkageCutoff):
  // convert R --> C++
  terms(Rcpp::as<std::vector<std::string>>(r_terms)),
  geneIDs(Rcpp::as<std::vector<std::string>>(r_geneIDs)),
  n_terms(int(terms.size())),
  
  // initialize data structures
  distMatrix(n_terms, terms),
  adjList(n_terms),
  clusList(terms),
  
  // initialize metrics
  dm(DistanceMetric(distanceMetric, distanceCutoff)),
  lm(LinkageMethod(linkageMethod, linkageCutoff, this->distFct()))
  { // checks: ensure vectors are of same size
    if (terms.size() != geneIDs.size())
      throw std::invalid_argument("input vectors (terms, geneIDs) must be the same size");
  };
  void computeDistances();
  void filterSeeds(); // informally denoting (node, neighbors) =: seed
  void mergeClusters();
  
  static constexpr double SAME_TERM_DISTANCE = -99;
  
  // reference to fast distance getting function to pass around
  std::function<double(int, int)> distFct() {
    return [this](int t1, int t2) {
      return distMatrix.getDistance(t1, t2);
    }; }
  
  Rcpp::NumericMatrix export_dm() const {return distMatrix.export_r();};
  Rcpp::DataFrame export_cl() const {return clusList.export_r();};
  
  
private:
  std::unordered_set<int>  filterSeed(int node, std::unordered_set<int> neighbors);
  ClusterList::ClusterIt findBestMergePartner(
      ClusterList::ClusterIt it1, std::list<std::unordered_set<int>>& clusters
  );
  
  // essential variables
  std::vector<std::string> terms;
  std::vector<std::string> geneIDs;
  int n_terms;
  
  // data structures
  DistanceMatrix distMatrix;
  AdjacencyList adjList;
  ClusterList clusList;
  
  // metrics
  DistanceMetric dm;
  LinkageMethod lm;
};

#endif /* richCluster_h */
