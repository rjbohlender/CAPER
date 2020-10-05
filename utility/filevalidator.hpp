//
// Created by Bohlender,Ryan James on 8/25/20.
//

#ifndef PERMUTE_ASSOCIATE_FILEVALIDATOR_HPP
#define PERMUTE_ASSOCIATE_FILEVALIDATOR_HPP

#include "split.hpp"
#include <set>

class FileValidator {
  RJBUtil::Splitter<std::string> matrix_header;
  static const std::set<std::string> matrix_variant_types;

public:
	FileValidator() = default;
	void validate_matrix_line(RJBUtil::Splitter<std::string> &line, int lineno) const;
	void validate_cov_line(RJBUtil::Splitter<std::string> &line, int lineno) const;
  	void validate_ped_line(RJBUtil::Splitter<std::string> &line, int lineno) const;
	void validate_bed_line(RJBUtil::Splitter<std::string> &line, int lineno) const;
	void validate_weight_line(RJBUtil::Splitter<std::string> &line, int lineno) const;
	void set_matrix_header(const std::string &header);
	static std::string build_error_message(int lineno, const std::string &message);
};

#endif //PERMUTE_ASSOCIATE_FILEVALIDATOR_HPP
