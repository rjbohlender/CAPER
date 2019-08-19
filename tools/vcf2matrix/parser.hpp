//
// Created by Bohlender,Ryan James on 2019-07-31.
//

#ifndef PERMUTE_ASSOCIATE_PARSER_HPP
#define PERMUTE_ASSOCIATE_PARSER_HPP

#include "reference.hpp"
#include <string>

class Parser {
  enum class Type {
    SNV,
    Insertion,
    Deletion,
    MNP
  };

  struct Variant {
	std::string location;
    std::vector<int> data;
  };

  struct GeneContainer {
    std::string gene;
    std::vector<std::string> transcripts;
    std::map<std::string, std::vector<Variant>> data;
    long endpos;
  };

public:
  Parser(const std::string &path, const std::string &outpath, Reference ref);



private:
  Reference ref_;
  std::vector<GeneContainer> genes_; // Genes, transcripts, variants

  void parse(std::istream &is, std::ostream &os);
  void cull(const std::string &chr, long pos, std::ostream &os); // Output genes we've finished
  static Type variant_classifier(const std::string &ref, const std::string &alt);
  static int variant_converter(const std::string &v);
};

#endif //PERMUTE_ASSOCIATE_PARSER_HPP
