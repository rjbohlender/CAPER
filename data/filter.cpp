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
  int lineno = -1;
  std::vector<std::string> methods;

  while(std::getline(ifs, line)) {
    lineno++;
    RJBUtil::Splitter<std::string> splitter(line, ",");
    if (lineno == 0) {
      std::copy(splitter.begin() + 1, splitter.end(), std::back_inserter(methods));
      continue;
    }
    if (boost::starts_with(line, "#")) { // Skip commented lines
      continue;
    }
    RJBUtil::Splitter<std::string> variant(splitter[0], ":");
    for (int i = 1; i < splitter.size(); i++) {
      if (splitter[i] == "1") {
		method_type_map[methods[i]].insert(variant[1]);
	  }
	}
  }
}

bool Filter::allow_variant(const std::string &method, RJBUtil::Splitter<std::string> &variant) {
  return method_type_map[method].find(variant[static_cast<int>(Indices::type)]) != method_type_map[method].end() &&
  		 method_type_map[method].find(variant[static_cast<int>(Indices::function)]) != method_type_map[method].end();
}
