//
// Created by Bohlender,Ryan James on 2019-07-31.
//

#include <iostream>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <sstream>
#include "parser.hpp"
#include "../../utility/filesystem.hpp"
#include "../../utility/split.hpp"

Parser::Parser(const std::string &path, const std::string &outpath, Reference ref)
	: ref_(std::move(ref)) {
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
  std::ofstream os(outpath);
  if (!os.good()) {
    throw(std::runtime_error("Failed to open output file for writing."));
  }

  parse(is, os);
}

/*!
 * @brief
 * @param is
 *
 * Read each variant, check if it is part of a gene, assemble a list of variants in a given gene based on it's bounds,
 * and output the gene when the sorted VCF is read beyond the end of the gene.
 */
void Parser::parse(std::istream &is, std::ostream &os) { // Parse VCF
  std::string line;
  std::string chr;

  std::vector<std::string> header;
  header.push_back("Gene");
  header.push_back("Transcript");
  header.push_back("Location");

  while (std::getline(is, line)) {
	if (boost::starts_with(line, "##")) {
	  continue;
	} else if (boost::starts_with(line, "#C")) {
	  // Handle samples
	  RJBUtil::Splitter<std::string> splitter(line, " \t");
	  std::copy(splitter.begin() + 9, splitter.end(), std::back_inserter(header));

	  // Output header
	  os << header[0];
	  for (int i = 1; i < header.size(); i++) {
	    os << "\t" << header[i];
	  }
	  os << std::endl;
	  continue;
	} // Fallthrough
        RJBUtil::Splitter<std::string> splitter(line, " \t");
        chr = splitter.str(0);
        long pos = std::stol(splitter.str(1));

	if (!boost::starts_with(chr, "chr")) {
	  chr = "chr" + chr;
	}

	cull(chr, pos, os); // Output completed genes.

	Type type = variant_classifier(splitter[3], splitter[4]);

	std::stringstream location_stream;

	switch (type) {
	case Type::SNV: location_stream << chr << "-" << pos << "-" << pos << "-SNV";
	  break;
	case Type::Deletion: location_stream << chr << "-" << pos << "-" << pos << "-deletion";
	  break;
	case Type::Insertion: location_stream << chr << "-" << pos << "-" << pos << "-insertion";
	  break;
	case Type::MNP: location_stream << chr << "-" << pos << "-" << pos << "-complex_substitution";
	  break;
	}

	Variant var{
		location_stream.str()
	};
        std::transform(splitter.begin() + 9, splitter.end(), std::back_inserter(var.data),
                       [](auto v) {
                         return Parser::variant_converter(
                             RJBUtil::Splitter<std::string>::str(v));
                       });

	std::vector<std::shared_ptr<Gene>> match = ref_.get_gene(chr, pos);
	for (const auto &p : match) {
	  bool gene_exists = false;
	  for (auto &g : genes_) {
	    if (g.gene == p->gene_name) {
		  for (const auto &t : p->transcript) {
			if (p->positions[t].first <= pos && p->positions[t].second >= pos) {
			  g.data[t].push_back(var);
			}
		  }
		  gene_exists = true;
		}
	  }
	  if (!gene_exists) {
		// Find end position
		long max_pos = 0;
		for (const auto &m : p->positions) {
		  if (m.second.second > max_pos) {
		    max_pos = m.second.second;
		  }
		}
		// Construct new container
		GeneContainer gene {
	      p->gene_name,
	      p->transcript,
		  std::map<std::string, std::vector<Variant>>(),
	      max_pos
	    };
		// Add variant
		for (const auto &t : p->transcript) {
		  gene.data[t] = std::vector<Variant>();
		  if (p->positions[t].first <= pos && p->positions[t].second >= pos) {
		    gene.data[t].push_back(var);
		  }
		}
		// Add gene
		genes_.push_back(gene);
	  }
	}
  }

  cull(chr, HUGE_VAL, os); // Finalize
}

void Parser::cull(const std::string &chr, long pos, std::ostream &os) {
  std::vector<int> to_remove;
  int i = -1;
  for (auto &g : genes_) {
    i++;
    if (pos > g.endpos) {
	  to_remove.push_back(i);
      // Output because this gene is finished
      for (const auto &t : g.transcripts) {
        for (const auto &v : g.data[t]) {
          os << g.gene << "\t" << t << "\t" << v.location;
          for (const auto &gt : v.data) {
            os << "\t" << gt;
          }
          os << std::endl;
        }
      }
    }
  }

  // Remove genes in reverse order, do not invalidate iterators
  for(int i = to_remove.size() - 1; i >= 0; i--) {
    genes_.erase(genes_.begin() + to_remove[i]);
  }
}

Parser::Type Parser::variant_classifier(const std::string &ref, const std::string &alt) {
  if (alt == "-") {
	return Type::Deletion;
  } else if (ref.size() < alt.size()
	  && alt.find(",") == alt.size()) { // Alt is larger than reference and not multiallelic
	return Type::Insertion;
  } else if (ref == "-") {
	return Type::Insertion;
  } else if (ref.size() > alt.size()) {
	return Type::Deletion;
  } else if (ref.size() == 1 && alt.size() == 1 && alt != "-") {
	return Type::SNV;
  } else if (alt.size() > 1 && alt.find(",") != alt.size()) {
	return Type::MNP;
  }
  std::cerr << "Variant: " << ref << "\t" << alt << std::endl;
  throw (std::runtime_error("Unknown type."));
}

int Parser::variant_converter(const std::string &v) {
  if (v.size() < 3) {
	return 5; // Missing
  } else if (v.rfind(".", 3) != std::string::npos) {
	return 5;
  }
  // Handle MNPs by converting from bool
  int b1 = v[0] - '0';
  int b2 = v[2] - '0';
  return static_cast<int>(b1 > 0) + static_cast<int>(b2 > 0);
}

