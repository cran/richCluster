//
//  DistanceMatrix.cpp
//  richCluster
//
//  Created by Sarah on 6/2/25.
//

#include <stdio.h>
#include "DistanceMatrix.h"
#include <Rcpp.h>

std::size_t DistanceMatrix::getDistanceIndex(int t1, int t2) const {
  // DS-12: widen BEFORE multiplying.  (row_index * n_terms) in int arithmetic
  // overflows -- signed overflow, undefined behaviour -- once n_terms reaches
  // 46341, and a sanitizer build reports it there.  The result is identical
  // for every n that fits in memory.
  const std::size_t row_index = static_cast<std::size_t>(t1);
  const std::size_t col_index = static_cast<std::size_t>(t2);
  return (row_index * static_cast<std::size_t>(n_terms)) + col_index;
}

double DistanceMatrix::getDistance(int t1, int t2) const {
  const std::size_t i = DistanceMatrix::getDistanceIndex(t1, t2);
  return distances[i];
} 

void DistanceMatrix::setDistance(double distance, int t1, int t2) {
  const std::size_t i = DistanceMatrix::getDistanceIndex(t1, t2);
  distances[i] = distance;
}

// export utility to R
Rcpp::NumericMatrix DistanceMatrix::export_r() const {
  Rcpp::NumericMatrix dm(n_terms, n_terms);
  // unpack the vector
  //
  // SPEC-RC-007 / T1-07, OQ-1 Reading 1 (author ruling, 2026-08-25): the
  // exported diagonal is 1, the self-similarity every metric agrees on --
  // kappa(A, A), Jaccard(A, A) and Dice(A, A) are all 1 by definition.
  //
  // The change is made HERE, at export time, and NOT at the storage site.  The
  // in-memory matrix keeps richCluster::SAME_TERM_DISTANCE, which is what
  // tools/cpp-units/rc_cpp_units.cpp's StubMatrix reproduces so the linkage
  // functions under unit test see exactly what they see in production.  Storing
  // 1 instead would invalidate that stub, and test-cpp-units.R:24 forbids
  // editing the harness to make an item green.
  //
  // Nothing downstream reads the diagonal: computeDistances() `continue`s on
  // i == j before the cutoff test, so no self-edge enters the adjacency list,
  // and LinkageMethod skips *i == *j.  Off-diagonal entries are untouched.
  for (int i = 0; i < n_terms; ++i) {
    for (int j = 0; j < n_terms; ++j) {
      dm(i, j) = (i == j) ? 1.0 : getDistance(i, j);
    }
  }
  Rcpp::List dimnames = Rcpp::List::create(terms, terms);
  dm.attr("dimnames") = dimnames;
  return dm;
}
