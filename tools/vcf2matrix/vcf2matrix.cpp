//
// Created by Bohlender,Ryan James on 2019-07-30.
//

#include <iostream>

#include <boost/program_options.hpp>
#include "reference.hpp"
#include "parser.hpp"

namespace po = boost::program_options;

int main(int argc, char **argv) {

  po::options_description desc("Tool for conversion of VCF files to Matrix files.");
  po::variables_map vm;

  try {
	desc.add_options()
			("vcf,i",
			 po::value<std::string>()->required(),
			 "VCF file path.")
			("matrix,o",
			 po::value<std::string>()->required(),
			 "Matrix output path.")
			("ref,r",
			 po::value<std::string>()->required(),
			 "RefFlat file for gene annotation.")
			("help,h",
			 "Display usage message.");
	po::store(po::parse_command_line(argc, argv, desc), vm);
	if (vm.count("help")) {
	  std::cerr << desc << "\n";
	  return 1;
	}
	po::notify(vm);
  } catch (po::required_option &e) {
	std::cerr << "Missing required option:\n" << e.what() << "\n";
	std::cerr << desc << "\n";
	return 1;
  } catch (std::exception &e) {
	std::cerr << "Error: " << e.what() << "\n";
	return 1;
  }

  Reference ref(vm["ref"].as<std::string>());
  Parser(vm["vcf"].as<std::string>(), vm["matrix"].as<std::string>(), ref);
}