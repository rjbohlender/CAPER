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
  parse(is);
}

void Reference::parse(std::istream &ifs) {
  std::string line;

  while (std::getline(ifs, line)) {
	RJBUtil::Splitter<std::string> splitter(line, " \t");
	// Skip ill-formed lines.
	if (splitter.size() < 11) {
	  continue;
	}

	RJBUtil::Splitter<std::string> ex_start(splitter[exon_starti], ",");
	RJBUtil::Splitter<std::string> ex_end(splitter[exon_endi], ",");

	if (ex_start.size() != ex_end.size()) {
	  throw (std::runtime_error("Mismatched exon start and end positions."));
	}

	// If gene hasn't been read yet.
	if (gene_map_.count(splitter[genei]) <= 0) {
	  if (data_.count(splitter[chri]) == 0) {
		data_[splitter[chri]] = std::vector<std::shared_ptr<Gene>>();
	  }
	  data_[splitter[chri]].push_back(std::make_shared<Gene>());
	  gene_map_[splitter[genei]] = data_[splitter[chri]].back();

	  // Update fields
	  try {
		data_[splitter[chri]].back()->gene_name = splitter[genei];
		// Checking if a transcript already exists so as to avoid duplicates.
		if (std::find(data_[splitter[chri]].back()->transcript.begin(), data_[splitter[chri]].back()->transcript.end(), splitter[txi]) == data_[splitter[chri]].back()->transcript.end()) {
		  data_[splitter[chri]].back()->transcript.push_back(splitter[txi]);
		}
		data_[splitter[chri]].back()->chromosome = splitter[chri];
		data_[splitter[chri]].back()->strand = splitter[strandi];
		data_[splitter[chri]].back()->positions[splitter[txi]] =
			std::make_pair(std::stol(splitter[starti]), std::stol(splitter[endi]));
		data_[splitter[chri]].back()->cds[splitter[txi]] =
			std::make_pair(std::stol(splitter[cds_starti]), std::stol(splitter[cds_endi]));
		data_[splitter[chri]].back()->exons[splitter[txi]] = std::stol(splitter[exonsi]);
		data_[splitter[chri]].back()->exon_spans[splitter[txi]] = std::vector<std::pair<long, long>>();
		for (int i = 0; i < ex_start.size(); i++) {
		  data_[splitter[chri]].back()->exon_spans[splitter[txi]]
			  .push_back(std::make_pair(std::stol(ex_start[i]), std::stol(ex_end[i])));
		}
	  } catch(std::invalid_argument &e) {
	    std::cerr << e.what() << std::endl;
	    std::cerr << line << std::endl;
	    std::exit(1);
	  }
	} else { // Add transcript to existing gene
	  // Update fields
	  if (gene_map_[splitter[genei]]->gene_name != splitter[genei]) {
	    std::cerr << gene_map_[splitter[genei]]->gene_name << "\t" << splitter[genei] << std::endl;
	    std::cerr << data_[splitter[chri]].size() << std::endl;
		throw (std::runtime_error("Gene name mismatch."));
	  }
	  if (std::find(gene_map_[splitter[genei]]->transcript.begin(), gene_map_[splitter[genei]]->transcript.end(), splitter[txi]) == gene_map_[splitter[genei]]->transcript.end()) {
		gene_map_[splitter[genei]]->transcript.push_back(splitter[txi]);
	  }
	  if (gene_map_[splitter[genei]]->positions.count(splitter[txi]) == 0) {
		gene_map_[splitter[genei]]->positions[splitter[txi]] = std::make_pair(std::stol(splitter[starti]), std::stol(splitter[endi]));
	  }
	  if (gene_map_[splitter[genei]]->cds.count(splitter[txi]) == 0) {
		gene_map_[splitter[genei]]->cds[splitter[txi]] = std::make_pair(std::stol(splitter[cds_starti]), std::stol(splitter[cds_endi]));
	  }
	  if (gene_map_[splitter[genei]]->exons.count(splitter[txi]) == 0) {
		gene_map_[splitter[genei]]->exons[splitter[txi]] = std::stol(splitter[exonsi]);
	  }
	  if (gene_map_[splitter[genei]]->exon_spans.count(splitter[txi]) == 0) {
		for (int i = 0; i < ex_start.size(); i++) {
		  gene_map_[splitter[genei]]->exon_spans[splitter[txi]].push_back(std::make_pair(std::stol(ex_start[i]), std::stol(ex_end[i])));
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
