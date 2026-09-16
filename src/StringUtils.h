//
//  StringUtils.h
//  richCluster
//
//  Created by Sarah on 6/2/25.
//

#ifndef StringUtils_h
#define StringUtils_h

#include <stdio.h>

#include <string>
#include <sstream>
#include <vector>
#include <unordered_set>
#include <set>

class StringUtils {
public:
  // (fast) methods for splitting strings to vectors and unordered_sets (specifying delimiter)
  static std::vector<std::string> splitStringToVector(const std::string& input, const std::string& delimiter);
  static std::unordered_set<std::string> splitStringToUnorderedSet(const std::string& input, const std::string& delimiter);
   
  // (slow) overloaded methods for splitting strings using regex (no delimiter specified)
  static std::vector<std::string> splitStringToVector(const std::string& input);
  static std::unordered_set<std::string> splitStringToUnorderedSet(const std::string& input);
   
  // converting vectors and sets --> strings
  static std::string vectorToString(const std::vector<std::string>& vector, const std::string& delimiter);
  // C9: std::set<int> variant -- iterates ascending on every standard library
  static std::string setToString(const std::set<int>& set, const std::string& delimiter);
  // template method for converting any std::unordered_set<T> to string
  template <typename T>
  static std::string unorderedSetToString(const std::unordered_set<T>& set, const std::string& delimiter);
   
  // used for counting total # geneIDs
  //
  // SPEC-RC-007 / T0-06: the delimiter is a parameter, defaulted to ",".  The
  // default is load-bearing -- DavidClustering.cpp calls this with one argument
  // and is deliberately out of scope for T0-06 (spec.md SS-C.4), so it keeps
  // 1.0.2 behaviour without being edited.
  static int countUniqueElements(const std::vector<std::string>& stringifiedVector,
                                 const std::string& delimiter = ",");
   
};

#endif /* StringUtils_h */
