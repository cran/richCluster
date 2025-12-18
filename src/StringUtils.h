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
  // template method for converting any std::unordered_set<T> to string
  template <typename T>
  static std::string unorderedSetToString(const std::unordered_set<T>& set, const std::string& delimiter);
   
  // used for counting total # geneIDs
  static int countUniqueElements(const std::vector<std::string>& stringifiedVector);
   
};

#endif /* StringUtils_h */
