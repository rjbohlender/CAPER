//
// Created by Bohlender,Ryan James on 11/5/20.
//

#include <fstream>
#include <boost/algorithm/string/predicate.hpp>
#include "filter.hpp"
#include "../utility/indices.hpp"
#include "../utility/filesystem.hpp"

Filter::Filter(const std::string &file_path) {
  if (!check_file_exists(file_path)) {
	throw(std::runtime_error("ERROR: Path to whitelist is incorrect."));
  }
  std::string line;
  std::ifstream ifs(file_path);

  while(std::getline(ifs, line)) {
    RJBUtil::Splitter<std::string> splitter(line, ":");
    if (splitter.size() != 2) {
      throw(std::runtime_error("ERROR: Whitelist malformed. Expected TEST:type."));
    }
    if (boost::starts_with(line, "#")) { // Skip commented lines
      continue;
    }
    method_type_map[splitter[0]].insert(splitter[1]);
  }
}

bool Filter::allow_variant(const std::string &method, RJBUtil::Splitter<std::string> &variant) {
  return method_type_map[method].find(variant[static_cast<int>(Indices::type)]) != method_type_map[method].end() &&
  		 method_type_map[method].find(variant[static_cast<int>(Indices::function)]) != method_type_map[method].end();
}
