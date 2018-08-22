//
// Created by Bohlender,Ryan James on 8/14/18.
//

#include <fstream>

#include "filesystem.hpp"


bool check_file_exists(const std::string &path) {
  std::ifstream ifs(path);
  return ifs.good();
}
