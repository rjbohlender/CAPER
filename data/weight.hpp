//
// Created by Bohlender,Ryan James on 8/14/18.
//

#ifndef PERMUTE_ASSOCIATE_CASM_HPP
#define PERMUTE_ASSOCIATE_CASM_HPP

#include <string>
#include <unordered_map>
#include <sstream>


class Weight {
public:
  Weight() = default;
  explicit Weight(const std::string &ifile);
  explicit Weight(std::stringstream &ifile);

  [[nodiscard]] double get(const std::string &k) const;
  double get(const std::string &k);

  [[nodiscard]] bool empty() const;
  bool empty();

private:
    static constexpr int field_count = 9;
    static constexpr int chrom_index = 0;
    static constexpr int start_index = 1;
    static constexpr int end_index = 2;
    static constexpr int ref_index = 3;
    static constexpr int alt_index = 4;
    static constexpr int type_index = 5;
    static constexpr int gene_index = 6;
    static constexpr int transcript_index = 7;
    static constexpr int weight_index = 8;

  std::unordered_map<std::string, double> weights_;

};

#endif //PERMUTE_ASSOCIATE_CASM_HPP
