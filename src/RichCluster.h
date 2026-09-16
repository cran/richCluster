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
#include <set>
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
              std::string linkageMethod, double linkageCutoff,
              std::string geneDelim = ",",
              bool verbose = false):
  // convert R --> C++
  terms(Rcpp::as<std::vector<std::string>>(r_terms)),
  geneIDs(Rcpp::as<std::vector<std::string>>(r_geneIDs)),
  n_terms(int(terms.size())),
  // SPEC-RC-007 / T0-06.  Declared adjacent to n_terms below and initialised in
  // the SAME relative position here: C++ initialises members in DECLARATION
  // order, and a mismatch warns under -Wreorder.
  geneDelim(geneDelim),
  // Gate for the progress narration.  Declared and initialised in the same
  // relative position as geneDelim above, for the -Wreorder reason stated there.
  verbose(verbose),

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
  std::set<int>  filterSeed(int node, std::set<int> neighbors);   // C9: ordered sets
  ClusterList::ClusterIt findBestMergePartner(
      ClusterList::ClusterIt it1, std::list<std::set<int>>& clusters
  );
  
  // essential variables
  std::vector<std::string> terms;
  std::vector<std::string> geneIDs;
  int n_terms;
  std::string geneDelim;   // SPEC-RC-007 / T0-06
  bool verbose;            // false => the core runs silent (CRAN default)

  // data structures
  DistanceMatrix distMatrix;
  AdjacencyList adjList;
  ClusterList clusList;
  
  // metrics
  DistanceMetric dm;
  LinkageMethod lm;
};

#endif /* richCluster_h */
