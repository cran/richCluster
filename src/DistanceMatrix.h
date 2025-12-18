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
    distances.resize(n_terms * n_terms);
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
  int getDistanceIndex(int t1, int t2) const;
};

#endif /* DistanceMatrix_h */
