//
//  ClusterList.cpp
//  richCluster
//
//  Created by Sarah on 6/3/25.
//

#include <stdio.h>
#include <Rcpp.h>
#include "StringUtils.h"
#include "ClusterList.h"
#include <string>
#include <vector>

Rcpp::DataFrame ClusterList::export_r() const {
  std::vector<std::string> termIndicesColumn;
  std::vector<std::string> termNamesColumn;
  std::vector<int> clusterColumn;
  
  int n = 1;
  // iterate through ClusterList and convert to vectors
  for (const auto& clusterGroup : clusterList) {
    // turn the ints into a string
    std::string termIndicesString = StringUtils::unorderedSetToString(clusterGroup, ", ");
    termIndicesColumn.push_back(termIndicesString); // append to termIndices
     
    std::vector<std::string> clusterGroupTerms;
    for (const auto& term_index : clusterGroup) {
      std::string term = terms[term_index];
      clusterGroupTerms.push_back(term);
    } 
    // convert vector/unordered_set to one comma-delimited string
    std::string clusterGroupTerms_string = StringUtils::vectorToString(clusterGroupTerms, ", ");
    termNamesColumn.push_back(clusterGroupTerms_string); // append to termNames
     
    // append cluster number to clusterColumn
    // append cluster number to clusterColumn
    clusterColumn.push_back(n++);
  } 
  
  //cCreate and return a DataFrame using Rcpp
  return Rcpp::DataFrame::create(Rcpp::Named("Cluster") = clusterColumn,
                                 Rcpp::Named("TermNames") = termNamesColumn,
                                 Rcpp::Named("TermIndices") = termIndicesColumn);
}

void ClusterList::deduplicate() {
  std::unordered_set<std::string> seen;
  auto it = clusterList.begin();
  while (it != clusterList.end()) {
    // Create a canonical string representation
    std::vector<int> sorted(it->begin(), it->end());
    std::sort(sorted.begin(), sorted.end());
    std::string key;
    for (int id : sorted) key += std::to_string(id) + ",";
    
    // Check for duplicates
    if (seen.count(key)) {
      it = clusterList.erase(it);
    } else { 
      seen.insert(key);
      ++it;
    } 
  }
}
