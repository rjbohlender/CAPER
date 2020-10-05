//
// Created by Bohlender,Ryan James on 8/25/20.
//

#include "filevalidator.hpp"
#include <sstream>

const std::set<std::string> FileValidator::matrix_variant_types {
	"SNV",
	"insertion",
	"deletion",
	"SPDA",
	"complex_substitution"
};

void FileValidator::validate_matrix_line(RJBUtil::Splitter<std::string> &line, int lineno) const {
  if (line.size() < 3) {
	std::string msg = build_error_message(lineno, "ERROR: Matrix Line Validation -- Line appears to be truncated. Line not long enough.");
	throw(std::runtime_error(msg.c_str()));
  }
  RJBUtil::Splitter<std::string> loc_splitter(line[2], "-");
  if(loc_splitter.size() != 4) {
	std::string msg = build_error_message(lineno, "ERROR: Matrix Line Validation -- Malformed location. Format should be chromosome-start-end-type.");
    throw(std::runtime_error(msg.c_str()));
  }
  if(matrix_variant_types.find(loc_splitter[3]) == matrix_variant_types.end()) {
	std::string msg = build_error_message(lineno, "ERROR: Matrix Line Validation -- Variant type in location incorrect. Must be one of {SNV, insertion, deletion, SPDA, complex_substitution}.");
    throw(std::runtime_error(msg.c_str()));
  }
}

void FileValidator::validate_cov_line(RJBUtil::Splitter<std::string> &line, int lineno) const {
  if(line.size() < 2) {
	std::string msg = build_error_message(lineno, "ERROR: COV Line Validation -- Line appears to be truncated. Line not long enough.");
	throw(std::runtime_error(msg.c_str()));
  }
  for(auto it = line.begin() + 1; it != line.end(); it++) {
    try {
      std::stod(*it);
    } catch (std::exception &e) {
	  std::string msg = build_error_message(lineno, "ERROR: COV Line Validation -- Non-numeric value provided in covariates.");
	  throw(std::runtime_error(msg.c_str()));
    }
  }
}

void FileValidator::validate_ped_line(RJBUtil::Splitter<std::string> &line, int lineno) const {
  if(line.size() < 6) {
	std::string msg = build_error_message(lineno, "ERROR: PED Line Validation -- Line appears to be truncated. Line not long enough.");
	throw(std::runtime_error(msg.c_str()));
  }
}

void FileValidator::validate_bed_line(RJBUtil::Splitter<std::string> &line, int lineno) const {
  if(line.size() < 3) {
	std::string msg = build_error_message(lineno, "ERROR: BED Line Validation -- Line appears to be truncated. Line not long enough.");
	throw(std::runtime_error(msg.c_str()));
  }
}

void FileValidator::validate_weight_line(RJBUtil::Splitter<std::string> &line, int lineno) const {
  if(line.size() < 5) {
	std::string msg = build_error_message(lineno, "ERROR: Weight Line Validation -- Line appears to be truncated. Line not long enough.");
	throw(std::runtime_error(msg.c_str()));
  }
  if(matrix_variant_types.find(line[3]) == matrix_variant_types.end()) {
	std::string msg = build_error_message(lineno, "ERROR: Weight Line Validation -- Variant type in location incorrect. Must be one of {SNV, insertion, deletion, SPDA, complex_substitution}.");
	throw(std::runtime_error(msg.c_str()));
  }
  try {
	std::stod(line[4]);
  } catch (std::exception &e) {
	std::string msg = build_error_message(lineno, "ERROR: Weight Line Validation -- Non-numeric value provided in weights.");
	throw(std::runtime_error(msg.c_str()));
  }
}

void FileValidator::set_matrix_header(const std::string &header) {
  matrix_header = RJBUtil::Splitter<std::string>(header, "\t");
}

std::string FileValidator::build_error_message(int lineno, const std::string &message) {
  std::stringstream ss;
  ss << message;
  ss << " Line: " << lineno;
  std::string msg = ss.str();
  return msg;
}

