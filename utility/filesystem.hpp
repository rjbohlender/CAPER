//
// Created by Bohlender,Ryan James on 8/14/18.
//

#ifndef PERMUTE_ASSOCIATE_FILESYSTEM_HPP
#define PERMUTE_ASSOCIATE_FILESYSTEM_HPP

#include <boost/iostreams/filtering_streambuf.hpp>
#include <string>

bool check_file_exists(const std::string &path);
bool check_directory_exists(const std::string &path);
bool is_gzipped(const std::string &path);
bool is_zstd(const std::string &path);
bool make_directory(const std::string &path);
std::string get_executable_path();
std::string get_executable_directory();
std::string resolve_default_filter_whitelist_path(
    const std::string &executable_path = "");

#endif //PERMUTE_ASSOCIATE_FILESYSTEM_HPP
