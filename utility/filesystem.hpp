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

#endif //PERMUTE_ASSOCIATE_FILESYSTEM_HPP
