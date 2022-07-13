//
// Created by Bohlender,Ryan James on 8/14/18.
//

#include <boost/iostreams/filter/bzip2.hpp>
#include <fstream>
#include <sstream>
#include <sys/stat.h>

#include "filesystem.hpp"


bool check_file_exists(const std::string &path) {
  std::ifstream ifs(path);
  return ifs.good();
}

bool check_directory_exists(const std::string &path) {
  struct stat st;
  return (stat(path.c_str(), &st) == 0);
}

bool is_gzipped(const std::string &path) {
  char byte1;
  char byte2;
  std::ifstream isource(path, std::ios_base::binary);
  if(check_file_exists(path)) {
    isource.get(byte1);
    isource.get(byte2);

    return (byte1 == '\x1F' && byte2 == '\x8B');
  } else {
    throw(std::logic_error("File doesn't exist."));
  }
}
