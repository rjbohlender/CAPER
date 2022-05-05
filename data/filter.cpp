//
// Created by Bohlender,Ryan James on 11/5/20.
//

#include <fstream>
#include <boost/algorithm/string/predicate.hpp>
#include <sstream>
#include "filter.hpp"
#include "../utility/indices.hpp"
#include "../utility/filesystem.hpp"

Filter::Filter(const std::string &file_path) {
  if (!check_file_exists(file_path)) {
        std::stringstream ss;
        ss << "ERROR: Path to whitelist is incorrect. " << file_path << std::endl;
	throw(std::runtime_error(ss.str().c_str()));
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
      if (strcmp(splitter[i].c_str(), "1") == 0) {
		method_type_map[methods[i - 1]].insert(variant[1]);
	  }
	}
  }
}

bool Filter::allow_variant(const std::string &method, RJBUtil::Splitter<std::string> &variant) {
  return method_type_map[method].find(variant[static_cast<int>(Indices::type)]) != method_type_map[method].end() &&
  		 method_type_map[method].find(variant[static_cast<int>(Indices::function)]) != method_type_map[method].end();
}
