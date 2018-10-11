//
// Created by Bohlender,Ryan James on 8/10/18.
//

#define CATCH_CONFIG_MAIN

#include <catch2/catch.hpp>
#include <sstream>
#include <string>
#include <map>
#include <armadillo>

#include "../data/covariates.hpp"
#include "../data/gene.hpp"
#include "../statistics/methods.hpp"
#include "../data/bed.hpp"


TEST_CASE( "Data Construction & Methods" ) {
  std::stringstream test_data;
  std::stringstream test_cov;
  std::stringstream test_casm;

  // Header
  test_data << "Gene\tTranscript\tPosition\tcase1\tcase2\tcontrol1\tcontrol2\n";
  // First Transcript
  test_data << "test_gene\ttest_transcript1\tchr1-1-1-SNV\t0\t1\t0\t1\n";
  test_data << "test_gene\ttest_transcript1\tchr1-25-25-SNV\t1\t1\t1\t1\n";
  test_data << "test_gene\ttest_transcript1\tchr1-27-38-insertion\t0\t2\t2\t2\n";
  test_data << "test_gene\ttest_transcript1\tchr1-55-55-SNV\t2\t0\t0\t1\n";
  // Second Transcript
  test_data << "test_gene\ttest_transcript2\tchr1-1-1-SNV\t0\t1\t0\t1\n";
  test_data << "test_gene\ttest_transcript2\tchr1-25-25-SNV\t1\t1\t9\t1\n";
  test_data << "test_gene\ttest_transcript2\tchr1-27-38-insertion\t0\t2\t2\t2\n";
  test_data << "test_gene\ttest_transcript2\tchr1-55-55-SNV\t2\t0\t0\t1\n";

  // Phenotypes and Covariates
  test_cov << "1\t1\t0.5\n";
  test_cov << "1\t0.5\t0.5\n";
  test_cov << "0\t0.5\t1\n";
  test_cov << "0\t1.5\t1.5\n";

  test_casm << "chr1\t1\t1\tSNV\t0.5\n";
  test_casm << "chr1\t25\t25\tSNV\t1.5\n";
  test_casm << "chr1\t27\t38\tinsertion\t3.0\n";
  test_casm << "chr1\t55\t55\tSNV\t0.5\n";

  Covariates cov(test_cov);
  // Don't log values
  Weight casm(test_casm, false);

  // Variant Counts
  std::map<std::string, unsigned long> nvariants { {"test_transcript1", 4},
                                                   {"test_transcript2", 4} };

  Gene gene(test_data, cov.get_nsamples(), nvariants, casm);

  SECTION( "Data Construction" ) {

    arma::vec phenotypes {1, 1, 0, 0};
    // Includes a row for the intercept.
    arma::mat covariates { {1.0, 1.0, 1.0, 1.0},
                           {1.0, 0.5, 0.5, 1.5},
                           {0.5, 0.5, 1.0, 1.5} };

    arma::mat transcript1 { {0.0, 1.0, 2.0, 2.0},
                            {1.0, 1.0, 0.0, 0.0},
                            {0.0, 1.0, 0.0, 0.0}, // Swapped for minor allele
                            {1.0, 1.0, 0.0, 1.0} };
    arma::mat transcript2 { {0.0, 1.0, 2.0, 2.0},
                            {1.0, 1.0, 0.0, 0.0}, // Removed missing data
                            {0.0, 0.0, 0.0, 0.0}, // Swapped for minor allele
                            {1.0, 1.0, 0.0, 1.0} };

    arma::vec weights {0.5, 1.5, 3.0, 0.5};

    REQUIRE( cov.get_ncases() == 2 );
    REQUIRE( cov.get_nsamples() == 4 );
    REQUIRE( arma::all(cov.get_phenotype_vector() == phenotypes) );
    REQUIRE( arma::all(arma::all(cov.get_covariate_matrix() == covariates)) );

    REQUIRE( gene.get_gene() == "test_gene" );
    REQUIRE( gene.get_transcripts()[0] == "test_transcript1" );
    REQUIRE( gene.get_transcripts()[1] == "test_transcript2" );
    REQUIRE( gene.get_transcripts().size() == 2 );

    REQUIRE( arma::all(arma::all(gene.get_matrix("test_transcript1") == transcript1)) );
    REQUIRE( arma::all(arma::all(gene.get_matrix("test_transcript2") == transcript2)) );
    REQUIRE( arma::all(gene.get_weights("test_transcript1") == weights) );
  }

  SECTION( "BED" ) {
    std::stringstream test_bed;

    test_bed << "chr1\t1234\t1234\n";
    test_bed << "chr1\t2345\t2346\n";
    test_bed << "chr1\t1234\t1235\tAF\n";
    test_bed << "chr1\t1234\t1234\n";
    test_bed << "chr6\t36867370\t36867370\n";
    test_bed << "chr6\t36867370\t36867370\tAF\n";
    test_bed << "chr6\t36867370\t36867371\n";

    Bed bed(test_bed);

    REQUIRE( !bed.empty() );
    REQUIRE( bed.size() == 2 );  // Two chromosomes
    REQUIRE( bed.chromosome_count("chr1") == 2 ); // Only two ranges for chromosome 1
    REQUIRE( bed.check_variant("chr1", 1234) );
    REQUIRE( bed.check_variant("chr1", 1235) );
    REQUIRE( bed.check_variant("chr6", "36867370") );
    REQUIRE( bed.check_variant("chr6", "36867371") );
    REQUIRE_THROWS( bed.check_variant("chr7", "36867371") );
  }

  SECTION( "RANK FUNCTION" ) {
    arma::vec test_vec {3, 5, 5, 5, 5, 8};
    arma::vec correct {1, 3.5, 3.5, 3.5, 3.5, 6};

    REQUIRE( arma::all(rank(test_vec, "ascend") == correct) );
    REQUIRE( arma::all(rank(test_vec, "descend") == arma::reverse(correct)) );

    test_vec = {3, 5, 5, 5, 7, 8};
    correct = {1, 3, 3, 3, 5, 6};

    REQUIRE( arma::all(rank(test_vec, "descend") == arma::reverse(correct)) );
    REQUIRE( arma::all(rank(test_vec, "ascend") == correct) );
  }

  SECTION("CALPHA") {
    Methods methods("CALPHA");

    double calpha = methods.call("test_transcript1", gene, cov);
    REQUIRE(calpha == Approx(-2.000000));

    calpha = methods.call("test_transcript2", gene, cov);
    REQUIRE(calpha == Approx(-1.500000));
  }

  SECTION("CMC") {
    Methods methods("CMC");

    double cmc = methods.call("test_transcript1", gene, cov);
    REQUIRE(cmc == Approx(0.449827));

    cmc = methods.call("test_transcript2", gene, cov);
    REQUIRE(cmc == Approx(0.591716));
  }

  SECTION("SKAT") {
    Methods methods("SKAT");

    double skat = methods.call("test_transcript1", gene, cov);
    REQUIRE(skat == Approx(0.625000));

    skat = methods.call("test_transcript2", gene, cov);
    REQUIRE(skat == Approx(0.750000));
  }

  SECTION("VAAST") {
    arma::vec log_casm;
    arma::vec test_pheno {1, 1, 1, 0 ,0};
    arma::mat test_geno { {1, 0},
                          {2, 1},
                          {0, 0},
                          {0, 1},
                          {0, 0} };

    Methods methods("VAAST");
    double test_vaast = methods.call(test_geno, test_pheno, log_casm);
    REQUIRE(test_vaast == Approx(1.8995198663785278));
    log_casm.reset();

    test_geno = { {1, 0},
                  {1, 1},
                  {1, 0},
                  {1, 1},
                  {1, 0} };
    test_vaast = methods.call(test_geno, test_pheno, log_casm, false, 0, 0);
    REQUIRE(test_vaast == Approx(0.10263280741763481));
    log_casm.reset();

    test_geno = { {1, 1},
                  {1, 1},
                  {1, 1},
                  {1, 1},
                  {1, 1} };
    log_casm = { 3.0, 5.0 };
    test_vaast = methods.call(test_geno, test_pheno, log_casm, false, 0, 0);
    REQUIRE(test_vaast == Approx(16.0));
  }

  SECTION("VT") {
    Methods methods("VT");

    double vt = methods.call("test_transcript1", gene, cov);
    REQUIRE(vt == Approx(0.452267));

    vt = methods.call("test_transcript2", gene, cov);
    REQUIRE(vt == Approx(0.408248));
  }

  SECTION("WSS") {
    Methods methods("WSS");

    double wss = methods.call("test_transcript1", gene, cov);
    REQUIRE(wss == 6.00);

    wss = methods.call("test_transcript2", gene, cov);
    REQUIRE(wss == 6.00);
  }
}
