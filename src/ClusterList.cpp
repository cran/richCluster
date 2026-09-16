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
#include <algorithm>
#include <string>
#include <vector>

Rcpp::DataFrame ClusterList::export_r() const {
  std::vector<std::string> termIndicesColumn;
  std::vector<std::string> termNamesColumn;
  std::vector<int> clusterColumn;
  
  // C9 (converge ledger, 2026-09-04): emit in a CONTENT-DEFINED order.  The
  // std::list order is the merge sequence; ordering the export by the member
  // sets themselves (std::set::operator< is lexicographic over ascending
  // members, so "smallest member first", full set as the tie-break) makes the
  // cluster numbering a function of the result alone -- identical on every
  // platform and for every permutation of the input.
  std::vector<const Cluster*> ordered;
  ordered.reserve(clusterList.size());
  for (const auto& c : clusterList) ordered.push_back(&c);
  std::stable_sort(ordered.begin(), ordered.end(),
                   [](const Cluster* a, const Cluster* b) { return *a < *b; });

  int n = 1;
  for (const Cluster* cp : ordered) {
    const Cluster& clusterGroup = *cp;
    // members are already ascending (std::set) -- canonical term order
    std::string termIndicesString = StringUtils::setToString(clusterGroup, ", ");
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
