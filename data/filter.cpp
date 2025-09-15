//
// Created by Bohlender,Ryan James on 11/5/20.
//

#include "filter.hpp"
#include "../utility/filesystem.hpp"
#include "matrix_indices.hpp"
#include <boost/algorithm/string/predicate.hpp>
#include <fstream>
#include <sstream>

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

  while (std::getline(ifs, line)) {
    lineno++;
    RJBUtil::Splitter<std::string> splitter(line, ",");
    if (lineno == 0) {
      for (auto it = splitter.begin() + 1; it != splitter.end(); ++it) {
        methods.emplace_back(std::string(*it));
      }
      continue;
    }
    if (boost::starts_with(line, "#")) { // Skip commented lines
      continue;
    }
    RJBUtil::Splitter<std::string> variant(splitter[0], ":");
    for (int i = 1; i < splitter.size(); i++) {
      if (splitter[i] == "1") {
        method_type_map[methods[i - 1]].insert(std::string(variant[1]));
      }
    }
  }
}

bool Filter::allow_variant(const std::string &method,
                           RJBUtil::Splitter<std::string> &variant) {
  auto type = variant.str(static_cast<int>(Indices::type));
  auto func = variant.str(static_cast<int>(Indices::function));
  return method_type_map[method].find(type) != method_type_map[method].end() &&
         method_type_map[method].find(func) != method_type_map[method].end();
}
