//
// Created by Bohlender,Ryan James on 8/14/18.
//

#include <boost/iostreams/filter/bzip2.hpp>
#include <array>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <sys/stat.h>
#include <unistd.h>

#include "filesystem.hpp"
#include "split.hpp"

#ifdef __APPLE__
#include <mach-o/dyld.h>
#endif

namespace fs = std::filesystem;

namespace {
std::string path_to_string(const fs::path &path) {
  return path.lexically_normal().string();
}

std::string directory_with_separator(const fs::path &path) {
  auto normalized = path_to_string(path);
  if (!normalized.empty() && normalized.back() != fs::path::preferred_separator) {
    normalized.push_back(fs::path::preferred_separator);
  }
  return normalized;
}
} // namespace

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

std::string get_executable_path() {
  char pathbuf[4096] = {0};

#ifdef __APPLE__
  uint32_t pathbufsize = sizeof(pathbuf);
  if (_NSGetExecutablePath(pathbuf, &pathbufsize) != 0) {
    throw std::runtime_error("Unable to determine executable path.");
  }
#else
  const ssize_t len = readlink("/proc/self/exe", pathbuf, sizeof(pathbuf) - 1);
  if (len == -1) {
    throw std::runtime_error("Unable to determine executable path.");
  }
  pathbuf[len] = '\0';
#endif

  return path_to_string(fs::path(pathbuf));
}

std::string get_executable_directory() {
  return directory_with_separator(fs::path(get_executable_path()).parent_path());
}

std::string resolve_default_filter_whitelist_path(
    const std::string &executable_path) {
  const fs::path executable =
      executable_path.empty() ? fs::path(get_executable_path())
                              : fs::path(executable_path);
  const fs::path executable_directory = executable.parent_path();
  const std::array<fs::path, 4> candidates = {
      executable_directory.parent_path() / "filter" / "filter_whitelist.csv",
      executable_directory / "filter" / "filter_whitelist.csv",
      fs::path("/filter/filter_whitelist.csv"),
      fs::current_path() / "filter" / "filter_whitelist.csv",
  };

  for (const auto &candidate : candidates) {
    if (check_file_exists(candidate.string())) {
      return path_to_string(candidate);
    }
  }

  return path_to_string(candidates.front());
}
