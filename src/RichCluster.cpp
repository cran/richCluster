//
//  RichCluster.cpp
//  RichCluster
//
//  Created by Sarah on 6/2/25.
//

#include <stdio.h>
#include <string>
#include "RichCluster.h"
#include "StringUtils.h"
#include <Rcpp.h>

void richCluster::computeDistances() {
  Rcpp::Rcout << "Computing distances..." << std::endl;
  int totalGeneCount = StringUtils::countUniqueElements(geneIDs);
  
  for (int i=0; i<n_terms; ++i) {
    // the unordered set of term1 genes
    std::unordered_set<std::string> term1_genes = StringUtils::splitStringToUnorderedSet(geneIDs[i], ",");
    
    for (int j=0; j<n_terms; ++j) {
      if (i == j) {
        distMatrix.setDistance(richCluster::SAME_TERM_DISTANCE, i, j);
        continue;
      }
      // unordered set of term2 genes
      std::unordered_set<std::string> term2_genes = StringUtils::splitStringToUnorderedSet(geneIDs[j], ",");
      
      double distanceScore = dm.computeDistance(term1_genes, term2_genes, totalGeneCount);
      distMatrix.setDistance(distanceScore, i, j);
      
      // if term similarity is ABOVE the threshold
      if (distanceScore >= dm.getCutoff()) {
        // add to adjacency list bidirectionally
        adjList.addNeighbor(i, j);
        adjList.addNeighbor(j, i);
      }
    }
  }
  Rcpp::Rcout << "Done filling out DistanceMatrix." << std::endl;
}

// go through adjacency list and find the best subset of each seed
void richCluster::filterSeeds() {
  Rcpp::Rcout << "Filtering seeds..." << std::endl;
  
  for (const auto& [node, neighbors] : richCluster::adjList.getAdjList()) {
    std::unordered_set<int> neighbors_set{neighbors};
    std::unordered_set<int> cluster = filterSeed(node, neighbors_set);
    clusList.addCluster(cluster);
  }
  Rcpp::Rcout << "Done filtering." << std::endl;
}

std::unordered_set<int> richCluster::filterSeed(
    int node, std::unordered_set<int> neighbors
) {
  std::unordered_set<int> cluster{node};
  while (true) {
    int bestN = -1;
    double bestLink = -1.0;
    
    for (int n : neighbors) {
      if (cluster.count(n)) continue;
      std::unordered_set<int> n_set{n};
      
      try {
        double link = lm.computeLinkage(cluster, n_set);
        if (link > bestLink) {
          bestLink = link;
          bestN = n;
        } 
      } catch (const std::exception& e) {
        Rcpp::Rcerr << "  EXCEPTION during linkage: " << e.what() << std::endl;
        throw; // rethrow to bubble up
      } 
    }
    if (bestLink < lm.getCutoff() || bestN == -1)
      break;
    cluster.insert(bestN);
  }
  return cluster;
}


void richCluster::mergeClusters() {
  Rcpp::Rcout << "Starting cluster merging..." << std::endl;

  bool mergingPossible = true;
  int iteration = 0;

  while (mergingPossible) {
    iteration++;
    Rcpp::Rcout << "Merge iteration " << iteration << "..." << std::endl;
    int nMerged = 0;
    auto& clusters = clusList.getList();
    
    for (auto it1 = clusters.begin(); it1 != clusters.end(); ++it1) {
      auto it2 = findBestMergePartner(it1, clusters);
      
      if (it2 != clusters.end() && it2 != it1) {
        // Merge cluster2 into cluster1
        it1->insert(it2->begin(), it2->end());
        clusters.erase(it2);  // immediately erase
        nMerged++;
      }
    } 
    Rcpp::Rcout << "  Number of merges in this iteration: " << nMerged << std::endl;
    if (nMerged == 0) {
      Rcpp::Rcout << "No more merges possible. Merging complete." << std::endl;
      break;
    } 
  }
  clusList.deduplicate();
}

ClusterList::ClusterIt richCluster::findBestMergePartner(
    ClusterList::ClusterIt it1, std::list<std::unordered_set<int>>& clusters
) {
  double bestLink = -1.0;
  auto bestIt = clusters.end();
   
  for (auto it2 = clusters.begin(); it2 != clusters.end(); ++it2) {
    if (it1 == it2) continue;
     
    double link = lm.computeLinkage(*it1, *it2);
    if (link > bestLink && link > lm.getCutoff()) {
      bestLink = link;
      bestIt = it2;
    } 
  }
  
  return bestIt;
} 



// the exported function to R
// [[Rcpp::export]]
Rcpp::List runRichCluster(Rcpp::CharacterVector terms,
                          Rcpp::CharacterVector geneIDs,
                          std::string distanceMetric, double distanceCutoff,
                          std::string linkageMethod, double linkageCutoff) {
  Rcpp::Rcout << "Starting richCluster..." << std::endl;
  Rcpp::Rcout << "terms.size = " << terms.size() << std::endl;
  Rcpp::Rcout << "geneIDs.size = " << geneIDs.size() << std::endl;
  try {
    richCluster RC(terms, geneIDs,
                   distanceMetric, distanceCutoff,
                   linkageMethod, linkageCutoff);
    RC.computeDistances();
    RC.filterSeeds();
    RC.mergeClusters();
    
    return Rcpp::List::create(
      Rcpp::_["distance_matrix"] = RC.export_dm(),
      Rcpp::_["all_clusters"]    = RC.export_cl()
    ); 
  } catch (const std::exception& e) {
    Rcpp::stop("C++ exception: %s", e.what());
  } catch (...) { 
    Rcpp::stop("Unknown C++ exception occurred.");
  } 
}
