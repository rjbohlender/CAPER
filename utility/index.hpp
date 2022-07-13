//
// Created by Bohlender,Ryan James on 10/7/20.
//

#ifndef PERMUTE_ASSOCIATE_INDEX_HPP
#define PERMUTE_ASSOCIATE_INDEX_HPP

#include "../utility/filesystem.hpp"
#include "../utility/split.hpp"
#include <zstr.hpp>
#include <fstream>
#include <iostream>
#include <map>

using namespace RJBUtil;

class Index {
  std::map<std::string, long> gene_map;
public:
  Index() = default;
  explicit Index(const std::string& fpath) {
    zstr::ifstream ifs_(fpath);

    std::string line;
    long lineno = 0;
    while(std::getline(ifs_, line)) {
      lineno++;
      Splitter<std::string> splitter(line, "\t ");
      if (splitter.size() != 2) {
        std::cerr << "Line " << lineno << " is malformed in the index file.\n";
        std::cerr << line << std::endl;
        std::exit(1);
      }
      try {
        gene_map[splitter[0]] = std::stol(splitter[1]);
      } catch(std::invalid_argument &e) {
        std::cerr << "Couldn't convert file position in index to long.\n";
        std::cerr << line;
        std::exit(1);
      }
    }
  }

  bool empty() {
    return gene_map.empty();
  }

  long gene_lookup(const std::string& gene_symbol) {
    if(gene_map.find(gene_symbol) == gene_map.end()) {
      throw(std::invalid_argument("Index doesn't contain the gene.\n"));
    }
    return gene_map[gene_symbol];
  }
};

#endif // PERMUTE_ASSOCIATE_INDEX_HPP
