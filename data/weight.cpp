//
// Created by Bohlender,Ryan James on 8/21/18.
//

#include "weight.hpp"
#include "../utility/filesystem.hpp"
#include "../utility/filevalidator.hpp"
#include "../utility/split.hpp"

#include <boost/algorithm/string/predicate.hpp>
#include <exception>
#include <stdexcept>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

/**
 * @brief Constructs a Weight object from a file
 *
 * This constructor reads weight data from the specified file path.
 * Each line should be tab-separated and contain chromosome, position, reference,
 * alternate, type, gene, transcript, and weight information.
 * Lines starting with '#' or empty lines are skipped.
 * For SNP type variants, they are stored as SNV in the internal map.
 *
 * @param ifile Path to the file containing weight data
 * @throws std::exception If weight value cannot be converted to double
 */
Weight::Weight(const std::string &ifile) {
  using namespace RJBUtil;
  if (!check_file_exists(ifile)) {
    std::cerr << "No weights provided." << std::endl;
    return;
  }
  std::ifstream ifs(ifile);
  std::string line;
  int lineno = -1;

  while (std::getline(ifs, line)) {
    lineno++;
    if (line.empty() || boost::starts_with(line, "#")) {
      continue;
    }

    Splitter<std::string> splitter(line, "\t");
    FileValidator::validate_weight_line(splitter, lineno);

    std::stringstream ss;
    // Example line:
    // chr13   114326148       114326148       G       A       SNV     CHAMP  NM_032436       2.55051363
    if (splitter[type_index] == "SNP") {
      ss << splitter[chrom_index] << ","
         << splitter[start_index] << ","
         << splitter[end_index] << ","
         << splitter[ref_index] << ","
         << splitter[alt_index] << ",SNV,"
         << splitter[gene_index] << ","
         << splitter[transcript_index];
    } else {
      ss << splitter[chrom_index] << ","
         << splitter[start_index] << ","
         << splitter[end_index] << ","
         << splitter[ref_index] << ","
         << splitter[alt_index] << ","
         << splitter[type_index] << ","
         << splitter[gene_index] << ","
         << splitter[transcript_index];
    }

    double weight = 1;
    try {
      weight = std::stod(splitter.at_str(weight_index));
    } catch (std::exception &e) {
      std::cerr << "Failed to convert weight to double. Line was: " << line
                << std::endl;
      std::cerr << "Line should be tab separated and formatted as <chrom> "
                   "<start_pos> <end_pos> <ref> <alt> <type> <gene> <transcript> <weight>"
                << std::endl;
      throw(e);
    }

    // Prevent math errors
    weights_[ss.str()] = weight;
  }
}

double Weight::get(const std::string &k) const {
  try {
    return weights_.at(k);
  } catch (std::exception &e) {
    std::cerr << "Failed to find weight: " << k << std::endl;
    return 1;
  }
}

double Weight::get(const std::string &k) { return weights_.at(k); }

bool Weight::empty() const { return weights_.empty(); }

bool Weight::empty() { return weights_.empty(); }

/**
 * @brief Constructs a Weight object from a stringstream
 *
 * This constructor reads weight data line by line from the provided stringstream.
 * Each line should be tab-separated and contain chromosome, position, reference,
 * alternate, type, gene, transcript, and weight information.
 * Lines starting with '#' or empty lines are skipped.
 * For SNP type variants, they are stored as SNV in the internal map.
 *
 * @param ifile stringstream containing the weight data in tab-separated format
 * @throws std::runtime_error If a line has incorrect format or field count
 * @throws std::exception If weight value cannot be converted to double
 */
Weight::Weight(std::stringstream &ifile) {
  std::string line;

  while (std::getline(ifile, line, '\n')) {
    // Skip empty lines and comments
    if (line.empty() || line[0] == '#') {
      continue;
    }
    RJBUtil::Splitter<std::string> splitter(line, "\t");

    // Validate line format
    if (splitter.size() != field_count) {
      std::cerr << "Incorrectly formatted weight line. Line was: " << line
                << std::endl;
      std::cerr << "Line should be tab separated and formatted as <chrom> "
                   "<start_pos> <end_pos> <ref> <alt> <type> <gene> <transcript> <weight>"
                << std::endl;
      throw(std::runtime_error("Incorrect line in weight file."));
    }

    // Build the key for the weights map
    std::stringstream ss;

    // Special handling for SNP type - convert to SNV
    if (splitter[type_index] == "SNP") {
      ss << splitter[chrom_index] << ","
         << splitter[start_index] << ","
         << splitter[end_index] << ","
         << splitter[ref_index] << ","
         << splitter[alt_index] << ",SNV,"
         << splitter[gene_index] << ","
         << splitter[transcript_index];
    } else {
      ss << splitter[chrom_index] << ","
         << splitter[start_index] << ","
         << splitter[end_index] << ","
         << splitter[ref_index] << ","
         << splitter[alt_index] << ","
         << splitter[type_index] << ","
         << splitter[gene_index] << ","
         << splitter[transcript_index];
    }

    // Parse and store the weight value
    double weight = 1;
    try {
      weight = std::stod(splitter.at_str(weight_index));
    } catch (std::exception &e) {
      std::cerr << "Failed to convert weight to double. Line was: " << line
                << std::endl;
      std::cerr << "Line should be tab separated and formatted as <chrom> "
                   "<start_pos> <end_pos> <ref> <alt> <type> <gene> <transcript> <weight>"
                << std::endl;
      throw(e);
    }

    // Store the weight in the map
    weights_[ss.str()] = weight;
  }
}
