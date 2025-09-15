//
// Created by Bohlender,Ryan James on 8/14/18.
//

#include <boost/iostreams/filter/bzip2.hpp>
#include <fstream>
#include <sstream>
#include <sys/stat.h>

#include "filesystem.hpp"
#include "split.hpp"

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

bool is_zstd(const std::string &path) {
  const uint8_t zstdref[4] = {0x28, 0xB5, 0x2F, 0xFD};
  uint8_t magic[4];
  std::ifstream isource(path, std::ios_base::binary);
  if(check_file_exists(path)) {
    isource.read((char *)magic, sizeof(magic));

    return (memcmp(magic, zstdref, sizeof(zstdref)) == 0);
  } else {
    throw(std::logic_error("File doesn't exist."));
  }
}

bool make_directory(const std::string &path) {
  // Split path into parts
  RJBUtil::Splitter<std::string> parts(path, "/");
  std::stringstream path_builder;
  // Absolute path
  if (path[0] == '/') {
    path_builder << "/";
  } else { // Relative path
    path_builder << "./";
  }
  for (const auto &part : parts) {
    path_builder << part << "/";
    if (!check_directory_exists(path_builder.str())) {
      // Read Write Exectute mask S_IRWXU
      int rc = mkdir(path_builder.str().c_str(), S_IRWXU);
      if (rc) {
        return false;
      }
    }
  }
  return true;
}
