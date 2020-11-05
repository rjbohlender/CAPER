//
// Created by Bohlender,Ryan James on 11/5/20.
//

#ifndef PERMUTE_ASSOCIATE_FILTER_HPP
#define PERMUTE_ASSOCIATE_FILTER_HPP

#include <string>
#include <map>
#include <set>
#include "../utility/split.hpp"

class Filter {
  std::map<std::string, std::set<std::string>> method_type_map;
public:
  explicit Filter(const std::string &file_path);

  bool allow_variant(const std::string &method, RJBUtil::Splitter<std::string> &variant);

};

#endif //PERMUTE_ASSOCIATE_FILTER_HPP
