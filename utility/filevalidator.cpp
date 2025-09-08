//
// Created by Bohlender,Ryan James on 8/25/20.
//

#include "filevalidator.hpp"
#include "../data/matrix_indices.hpp"
#include "split.hpp"
#include <cstddef>
#include <exception>
#include <sstream>
#include <string>
#include <stdexcept>
#include <iostream>
#include <set>

const std::set<std::string> FileValidator::matrix_variant_types{
    "SNV", "insertion", "deletion", "SPDA", "complex_substitution"};

void FileValidator::validate_matrix_line(RJBUtil::Splitter<std::string> &line,
                                         int lineno) const {
  if (line.size() < static_cast<int>(Indices::first)) {
    const std::string msg =
        build_error_message("ERROR: Matrix Line Validation -- Line appears to "
                            "be truncated. Line not long enough.",
                            lineno);
    throw(std::runtime_error(msg.c_str()));
  }
  if (!matrix_variant_types.contains(line[static_cast<int>(Indices::type)])
  ) {
    const std::string msg = build_error_message(
        "ERROR: Matrix Line Validation -- Variant type incorrect. Must be one "
        "of {SNV, insertion, deletion, SPDA, complex_substitution}.",
        lineno);
    throw(std::runtime_error(msg.c_str()));
  }
  if (line[static_cast<int>(Indices::first)].size() != matrix_sample_count) {
    const std::string msg =
        build_error_message("ERROR: Matrix Line Validation -- Genotype size "
                            "differs from sample count.",
                            lineno);
    throw(std::runtime_error(msg.c_str()));
  }
}

void FileValidator::validate_cov_line(RJBUtil::Splitter<std::string> &line,
                                      int lineno) {
  if (line.size() < cov_line_size) {
    const std::string msg =
        build_error_message("ERROR: COV Line Validation -- Line appears to be "
                            "truncated. Line not long enough."
                            "\nExpected format: sample_ID\tcovariate1\t...",
                            lineno);
    throw(std::runtime_error(msg.c_str()));
  }
}

void FileValidator::validate_ped_line(RJBUtil::Splitter<std::string> &line,
                                      int lineno) {
  if (line.size() < ped_line_size) {
    const std::string msg =
        build_error_message("ERROR: PED Line Validation -- Line appears to be "
                            "truncated. Line not long enough.",
                            lineno);
    throw(std::runtime_error(msg.c_str()));
  }
}

void FileValidator::validate_bed_line(RJBUtil::Splitter<std::string> &line,
                                      int lineno) {
  if (line.size() < bed_line_size) {
    const std::string msg = build_error_message(
        "ERROR: BED Line Validation -- Line appears to be truncated. Line not "
        "long enough."
        "\nExpected format: chrom\tstart_pos\tend_pos\tref\talt\t...",
        lineno);
    throw(std::runtime_error(msg.c_str()));
  }
}

void FileValidator::validate_weight_line(RJBUtil::Splitter<std::string> &line,
                                         int lineno) {
  if (line.size() != weight_line_size) {
    const std::string msg =
        build_error_message("ERROR: Weight Line Validation -- Line appears to "
                            "be truncated. Line not long enough."
                            "Line should be tab separated and formatted as <chrom> "
                            "<start_pos> <end_pos> <ref> <alt> <type> <gene> <transcript> <weight>",
                            lineno);
    throw(std::runtime_error(msg.c_str()));
  }
  if (!matrix_variant_types.contains(line[weight_type_index])) {
    const std::string msg = build_error_message(
        "ERROR: Weight Line Validation -- Variant type "
        "incorrect. Must be one of {SNV, insertion, deletion, "
        "SPDA, complex_substitution}.",
        lineno);
    throw(std::runtime_error(msg.c_str()));
  }
  try {
    std::stod(line.back());
  } catch (std::exception &e) {
    const std::string msg =
        build_error_message("ERROR: Weight Line Validation -- Non-numeric "
                            "value provided in weights.",
                            lineno);
    throw(std::runtime_error(msg.c_str()));
  }
}

void FileValidator::set_matrix_header(const std::string &header) {
  // Pass an owning string to Splitter so that the header fields referenced by
  // matrix_header remain valid even after the original `header` argument goes
  // out of scope. Using `std::string(header)` creates an rvalue that invokes
  // the owning constructor of Splitter.
  matrix_header =
      RJBUtil::Splitter<std::string>(std::string(header), "\t");
  matrix_sample_count =
      matrix_header.size() - static_cast<size_t>(Indices::first);
}

std::string FileValidator::build_error_message(const std::string &message,
                                               int lineno) {
  std::stringstream ss;
  ss << message;
  ss << " Line: " << lineno;
  std::string msg = ss.str();
  return msg;
}
