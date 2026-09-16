//
//  DistanceMatrix.h
//  richCluster
//
//  Created by Sarah on 6/2/25.
//

#ifndef DistanceMatrix_h
#define DistanceMatrix_h

#include <Rcpp.h>

class DistanceMatrix {
public:
  DistanceMatrix(int n_terms, std::vector<std::string>& terms):
  n_terms(n_terms), terms(terms) {
    // DS-12: n_terms * n_terms in int arithmetic is signed overflow -- and so
    // undefined behaviour -- from n_terms = 46341 upward.  Widen to size_t
    // before multiplying.  The cast changes nothing for any n that fits in
    // memory; it removes the UB that a sanitizer build flags.
    distances.resize(static_cast<std::size_t>(n_terms) *
                     static_cast<std::size_t>(n_terms));
  };
  
  double getDistance(int t1, int t2) const;
  void setDistance(double distance, int t1, int t2);
  
  Rcpp::NumericMatrix export_r() const;
  
private:
  std::vector<double> distances; // matrix is internally stored flattened
  
  // useful vars
  int n_terms;
  std::vector<std::string> terms;
  
  // index into flattened list
  std::size_t getDistanceIndex(int t1, int t2) const;
};

#endif /* DistanceMatrix_h */
