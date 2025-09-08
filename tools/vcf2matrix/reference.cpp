//
// Created by Bohlender,Ryan James on 2019-07-30.
//

#include <boost/iostreams/filter/gzip.hpp>
#include <iostream>
#include <algorithm>
#include "reference.hpp"

#include "../../utility/filesystem.hpp"
#include "../../utility/split.hpp"

Reference::Reference(const std::string &path) {
  std::ifstream ifs;
  boost::iostreams::filtering_streambuf<boost::iostreams::input> streambuf;

  if (is_gzipped(path)) {
	ifs.open(path, std::ios_base::in | std::ios_base::binary);
	streambuf.push(boost::iostreams::gzip_decompressor());
	streambuf.push(ifs);
  } else {
	ifs.open(path, std::ios_base::in);
	streambuf.push(ifs);
  }

  std::istream is(&streambuf);
  parse_refFlat(is);
}

void Reference::parse_refFlat(std::istream &ifs) {
  std::string line;

  while (std::getline(ifs, line)) {
	RJBUtil::Splitter<std::string> splitter(line, " \t");
	// Skip ill-formed lines.
	if (splitter.size() < 11) {
	  continue;
	}

	RJBUtil::Splitter<std::string> ex_start(splitter[static_cast<int>(refFlat::exon_starti)], ",");
	RJBUtil::Splitter<std::string> ex_end(splitter[static_cast<int>(refFlat::exon_endi)], ",");

	if (ex_start.size() != ex_end.size()) {
	  throw (std::runtime_error("Mismatched exon start and end positions."));
	}

	// If gene hasn't been read yet.
	if (gene_map_.count(splitter[static_cast<int>(refFlat::genei)]) <= 0) {
	  if (data_.count(splitter[static_cast<int>(refFlat::chri)]) == 0) {
		data_[splitter[static_cast<int>(refFlat::chri)]] = std::vector<std::shared_ptr<Gene>>();
	  }
	  data_[splitter[static_cast<int>(refFlat::chri)]].push_back(std::make_shared<Gene>());
	  gene_map_[splitter[static_cast<int>(refFlat::genei)]] = data_[splitter[static_cast<int>(refFlat::chri)]].back();

	  // Update fields
	  try {
		data_[splitter[static_cast<int>(refFlat::chri)]].back()->gene_name = splitter[static_cast<int>(refFlat::genei)];
		// Checking if a transcript already exists so as to avoid duplicates.
		if (std::find(data_[splitter[static_cast<int>(refFlat::chri)]].back()->transcript.begin(), data_[splitter[static_cast<int>(refFlat::chri)]].back()->transcript.end(), splitter[static_cast<int>(refFlat::txi)]) == data_[splitter[static_cast<int>(refFlat::chri)]].back()->transcript.end()) {
		  data_[splitter[static_cast<int>(refFlat::chri)]].back()->transcript.push_back(splitter[static_cast<int>(refFlat::txi)]);
		}
		data_[splitter[static_cast<int>(refFlat::chri)]].back()->chromosome = splitter[static_cast<int>(refFlat::chri)];
		data_[splitter[static_cast<int>(refFlat::chri)]].back()->strand = splitter[static_cast<int>(refFlat::strandi)];
                  data_[splitter[static_cast<int>(refFlat::chri)]].back()->positions[splitter[static_cast<int>(refFlat::txi)]] =
                          std::make_pair(std::stol(splitter.str(static_cast<int>(refFlat::starti))),
                                         std::stol(splitter.str(static_cast<int>(refFlat::endi))));
                  data_[splitter[static_cast<int>(refFlat::chri)]].back()->cds[splitter[static_cast<int>(refFlat::txi)]] =
                          std::make_pair(std::stol(splitter.str(static_cast<int>(refFlat::cds_starti))),
                                         std::stol(splitter.str(static_cast<int>(refFlat::cds_endi))));
                  data_[splitter[static_cast<int>(refFlat::chri)]].back()->exons[splitter[static_cast<int>(refFlat::txi)]] =
                      std::stol(splitter.str(static_cast<int>(refFlat::exonsi)));
                  data_[splitter[static_cast<int>(refFlat::chri)]].back()->exon_spans[splitter[static_cast<int>(refFlat::txi)]] = std::vector<std::pair<long, long>>();
                  for (int i = 0; i < ex_start.size(); i++) {
                    data_[splitter[static_cast<int>(refFlat::chri)]].back()->exon_spans[splitter[static_cast<int>(refFlat::txi)]]
                            .push_back(std::make_pair(std::stol(ex_start.str(i)), std::stol(ex_end.str(i))));
		}
	  } catch(std::invalid_argument &e) {
	    std::cerr << e.what() << std::endl;
	    std::cerr << line << std::endl;
	    std::exit(1);
	  }
	} else { // Add transcript to existing gene
	  // Update fields
	  if (gene_map_[splitter[static_cast<int>(refFlat::genei)]]->gene_name != splitter[static_cast<int>(refFlat::genei)]) {
	    std::cerr << gene_map_[splitter[static_cast<int>(refFlat::genei)]]->gene_name << "\t" << splitter[static_cast<int>(refFlat::genei)] << std::endl;
	    std::cerr << data_[splitter[static_cast<int>(refFlat::chri)]].size() << std::endl;
		throw (std::runtime_error("Gene name mismatch."));
	  }
	  if (std::find(gene_map_[splitter[static_cast<int>(refFlat::genei)]]->transcript.begin(), gene_map_[splitter[static_cast<int>(refFlat::genei)]]->transcript.end(), splitter[static_cast<int>(refFlat::txi)]) == gene_map_[splitter[static_cast<int>(refFlat::genei)]]->transcript.end()) {
		gene_map_[splitter[static_cast<int>(refFlat::genei)]]->transcript.push_back(splitter[static_cast<int>(refFlat::txi)]);
	  }
	  if (gene_map_[splitter[static_cast<int>(refFlat::genei)]]->positions.count(splitter[static_cast<int>(refFlat::txi)]) == 0) {
                  gene_map_[splitter[static_cast<int>(refFlat::genei)]]->positions[splitter[static_cast<int>(refFlat::txi)]] = std::make_pair(std::stol(splitter.str(static_cast<int>(refFlat::starti))), std::stol(splitter.str(static_cast<int>(refFlat::endi))));
	  }
	  if (gene_map_[splitter[static_cast<int>(refFlat::genei)]]->cds.count(splitter[static_cast<int>(refFlat::txi)]) == 0) {
                  gene_map_[splitter[static_cast<int>(refFlat::genei)]]->cds[splitter[static_cast<int>(refFlat::txi)]] = std::make_pair(std::stol(splitter.str(static_cast<int>(refFlat::cds_starti))), std::stol(splitter.str(static_cast<int>(refFlat::cds_endi))));
	  }
	  if (gene_map_[splitter[static_cast<int>(refFlat::genei)]]->exons.count(splitter[static_cast<int>(refFlat::txi)]) == 0) {
                  gene_map_[splitter[static_cast<int>(refFlat::genei)]]->exons[splitter[static_cast<int>(refFlat::txi)]] = std::stol(splitter.str(static_cast<int>(refFlat::exonsi)));
	  }
	  if (gene_map_[splitter[static_cast<int>(refFlat::genei)]]->exon_spans.count(splitter[static_cast<int>(refFlat::txi)]) == 0) {
		for (int i = 0; i < ex_start.size(); i++) {
                    gene_map_[splitter[static_cast<int>(refFlat::genei)]]->exon_spans[splitter[static_cast<int>(refFlat::txi)]].push_back(std::make_pair(std::stol(ex_start.str(i)), std::stol(ex_end.str(i))));
		}
	  }
	}
  }
  // Sort data_
  for (const auto &c : chromosomes) {
	std::sort(data_[c].begin(), data_[c].end(), [](auto &lhs, auto &rhs){ return (*lhs) < (*rhs); });
  }
}

std::vector<std::shared_ptr<Gene>> Reference::get_gene(const std::string &chr, long pos) {
  std::vector<std::shared_ptr<Gene>> ret;
  for (const auto &gp : data_[chr]) {
    if (*gp == pos) {
      ret.push_back(gp);
    }
  }

  return ret;
}
