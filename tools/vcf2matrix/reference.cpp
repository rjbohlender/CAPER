//
// Created by Bohlender,Ryan James on 2019-07-30.
//

#include <boost/iostreams/filter/gzip.hpp>
#include <iostream>
#include <algorithm>
#include "reference.hpp"

#include "utility/filesystem.hpp"
#include "utility/split.hpp"

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

        std::string gene = splitter.str(static_cast<int>(refFlat::genei));
        std::string chr = splitter.str(static_cast<int>(refFlat::chri));
        std::string tx = splitter.str(static_cast<int>(refFlat::txi));
        std::string strand = splitter.str(static_cast<int>(refFlat::strandi));

        // If gene hasn't been read yet.
        if (gene_map_.count(gene) <= 0) {
          if (data_.count(chr) == 0) {
                data_[chr] = std::vector<std::shared_ptr<Gene>>();
          }
          data_[chr].push_back(std::make_shared<Gene>());
          gene_map_[gene] = data_[chr].back();

          // Update fields
          try {
                data_[chr].back()->gene_name = gene;
                // Checking if a transcript already exists so as to avoid duplicates.
                if (std::find(data_[chr].back()->transcript.begin(), data_[chr].back()->transcript.end(), tx) == data_[chr].back()->transcript.end()) {
                  data_[chr].back()->transcript.push_back(tx);
                }
                data_[chr].back()->chromosome = chr;
                data_[chr].back()->strand = strand;
                  data_[chr].back()->positions[tx] =
                          std::make_pair(std::stol(splitter.str(static_cast<int>(refFlat::starti))),
                                         std::stol(splitter.str(static_cast<int>(refFlat::endi))));
                  data_[chr].back()->cds[tx] =
                          std::make_pair(std::stol(splitter.str(static_cast<int>(refFlat::cds_starti))),
                                         std::stol(splitter.str(static_cast<int>(refFlat::cds_endi))));
                  data_[chr].back()->exons[tx] =
                      std::stol(splitter.str(static_cast<int>(refFlat::exonsi)));
                  data_[chr].back()->exon_spans[tx] = std::vector<std::pair<long, long>>();
                  for (int i = 0; i < ex_start.size(); i++) {
                    data_[chr].back()->exon_spans[tx]
                            .push_back(std::make_pair(std::stol(ex_start.str(i)), std::stol(ex_end.str(i))));
                }
          } catch(std::invalid_argument &e) {
            std::cerr << e.what() << std::endl;
            std::cerr << line << std::endl;
            std::exit(1);
          }
        } else { // Add transcript to existing gene
          // Update fields
          if (gene_map_[gene]->gene_name != gene) {
            std::cerr << gene_map_[gene]->gene_name << "\t" << gene << std::endl;
            std::cerr << data_[chr].size() << std::endl;
                throw (std::runtime_error("Gene name mismatch."));
          }
          if (std::find(gene_map_[gene]->transcript.begin(), gene_map_[gene]->transcript.end(), tx) == gene_map_[gene]->transcript.end()) {
                gene_map_[gene]->transcript.push_back(tx);
          }
          if (gene_map_[gene]->positions.count(tx) == 0) {
                  gene_map_[gene]->positions[tx] = std::make_pair(std::stol(splitter.str(static_cast<int>(refFlat::starti))), std::stol(splitter.str(static_cast<int>(refFlat::endi))));
          }
          if (gene_map_[gene]->cds.count(tx) == 0) {
                  gene_map_[gene]->cds[tx] = std::make_pair(std::stol(splitter.str(static_cast<int>(refFlat::cds_starti))), std::stol(splitter.str(static_cast<int>(refFlat::cds_endi))));
          }
          if (gene_map_[gene]->exons.count(tx) == 0) {
                  gene_map_[gene]->exons[tx] = std::stol(splitter.str(static_cast<int>(refFlat::exonsi)));
          }
          if (gene_map_[gene]->exon_spans.count(tx) == 0) {
                for (int i = 0; i < ex_start.size(); i++) {
                    gene_map_[gene]->exon_spans[tx].push_back(std::make_pair(std::stol(ex_start.str(i)), std::stol(ex_end.str(i))));
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
