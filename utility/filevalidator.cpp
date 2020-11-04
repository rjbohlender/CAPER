//
// Created by Bohlender,Ryan James on 8/25/20.
//

#include "filevalidator.hpp"
#include "indices.hpp"
#include <sstream>

const std::set<std::string> FileValidator::matrix_variant_types {
	"SNV",
	"insertion",
	"deletion",
	"SPDA",
	"complex_substitution"
};

void FileValidator::validate_matrix_line(RJBUtil::Splitter<std::string> &line, int lineno) {
  if (line.size() < static_cast<int>(Indices::first)) {
	std::string msg = build_error_message(
		"ERROR: Matrix Line Validation -- Line appears to be truncated. Line not long enough.",
		lineno);
	throw(std::runtime_error(msg.c_str()));
  }
  if(matrix_variant_types.find(line[static_cast<int>(Indices::type)]) == matrix_variant_types.end()) {
	std::string msg = build_error_message(
		"ERROR: Matrix Line Validation -- Variant type incorrect. Must be one of {SNV, insertion, deletion, SPDA, complex_substitution}.",
		lineno);
    throw(std::runtime_error(msg.c_str()));
  }
}

void FileValidator::validate_cov_line(RJBUtil::Splitter<std::string> &line, int lineno) {
  if(line.size() < 2) {
	std::string msg = build_error_message(
		"ERROR: COV Line Validation -- Line appears to be truncated. Line not long enough.",
		lineno);
	throw(std::runtime_error(msg.c_str()));
  }
  for(auto it = line.begin() + 1; it != line.end(); it++) {
    try {
      std::stod(*it);
    } catch (std::exception &e) {
	  std::string msg = build_error_message("ERROR: COV Line Validation -- Non-numeric value provided in covariates.",
											lineno);
	  throw(std::runtime_error(msg.c_str()));
    }
  }
}

void FileValidator::validate_ped_line(RJBUtil::Splitter<std::string> &line, int lineno) {
  if(line.size() < 6) {
	std::string msg = build_error_message(
		"ERROR: PED Line Validation -- Line appears to be truncated. Line not long enough.",
		lineno);
	throw(std::runtime_error(msg.c_str()));
  }
}

void FileValidator::validate_bed_line(RJBUtil::Splitter<std::string> &line, int lineno) {
  if(line.size() < 3) {
	std::string msg = build_error_message(
		"ERROR: BED Line Validation -- Line appears to be truncated. Line not long enough.",
		lineno);
	throw(std::runtime_error(msg.c_str()));
  }
}

void FileValidator::validate_weight_line(RJBUtil::Splitter<std::string> &line, int lineno) {
  if(line.size() < 5) {
	std::string msg = build_error_message(
		"ERROR: Weight Line Validation -- Line appears to be truncated. Line not long enough.",
		lineno);
	throw(std::runtime_error(msg.c_str()));
  }
  if(matrix_variant_types.find(line[3]) == matrix_variant_types.end()) {
	std::string msg = build_error_message(
		"ERROR: Weight Line Validation -- Variant type in location incorrect. Must be one of {SNV, insertion, deletion, SPDA, complex_substitution}.",
		lineno);
	throw(std::runtime_error(msg.c_str()));
  }
  try {
	std::stod(line[4]);
  } catch (std::exception &e) {
	std::string msg = build_error_message("ERROR: Weight Line Validation -- Non-numeric value provided in weights.",
										  lineno);
	throw(std::runtime_error(msg.c_str()));
  }
}

void FileValidator::set_matrix_header(const std::string &header) {
  matrix_header = RJBUtil::Splitter<std::string>(header, "\t");
}

std::string FileValidator::build_error_message(const std::string &message, int lineno) {
  std::stringstream ss;
  ss << message;
  ss << " Line: " << lineno;
  std::string msg = ss.str();
  return msg;
}

