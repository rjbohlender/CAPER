//
// Created by Bohlender,Ryan James on 8/10/18.
//

#define CATCH_CONFIG_MAIN
#define CATCH_CONFIG_NO_POSIX_SIGNALS

#include <armadillo>
#include <catch2/catch.hpp>
#include <cmath>
#include <limits>
#include <map>
#include <sstream>
#include <string>

#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions/fisher_f.hpp>

#include "../data/bed.hpp"
#include "../data/covariates.hpp"
#include "../data/gene.hpp"
#include "../link/binomial.hpp"
#include "../statistics/glm.hpp"
#include "../statistics/methods.hpp"
#include "../utility/math.hpp"

namespace {
double cmc_dense_reference(const arma::mat &X, const arma::vec &Y,
                           double rare_freq, bool hotellings,
                           arma::uword nperm) {
  if (hotellings) {
    const double N = static_cast<double>(Y.n_rows);
    const double nA = arma::accu(Y);
    const double nU = N - nA;

    arma::rowvec maf = arma::sum(X, 0) / (2.0 * N);
    const arma::uvec rare = arma::find(maf < rare_freq);
    const arma::uvec common = arma::find(maf >= rare_freq);

    arma::mat Xnew;
    if (rare.n_elem <= 1) {
      Xnew = X;
    } else {
      arma::mat collapse = arma::sum(X.cols(rare), 1);
      collapse.transform([](double value) { return value > 1.0 ? 1.0 : value; });

      if (common.n_elem > 0) {
        Xnew = arma::join_horiz(X.cols(common), collapse);
      } else {
        Xnew = collapse;
      }
    }

    Xnew -= 1.0;

    arma::mat cases = Xnew.rows(arma::find(Y == 1));
    arma::mat controls = Xnew.rows(arma::find(Y == 0));

    arma::rowvec case_mean = arma::mean(cases);
    arma::rowvec control_mean = arma::mean(controls);

    arma::mat cov =
        ((nA - 1.0) * arma::cov(cases) + (nU - 1.0) * arma::cov(controls)) /
        (N - 2.0);
    arma::mat inv;
    if (!arma::inv_sympd(inv, cov)) {
      arma::pinv(inv, cov);
    }

    const double ret = arma::as_scalar((case_mean - control_mean) * inv *
                                       (case_mean - control_mean).t() * nA *
                                       nU / N);
    const double p = static_cast<double>(case_mean.n_elem);
    double stat = ret * (nA + nU - 1.0 - p) / (p * (nA + nU - 2.0));
    if (stat < 0.0) {
      stat = 0.0;
    }

    if (nperm > 0) {
      return stat;
    }

    boost::math::fisher_f fisher_f(p, nA + nU - 1.0 - p);
    if (std::isnan(stat)) {
      return 1.0;
    }

    return boost::math::cdf(boost::math::complement(fisher_f, stat));
  }

  arma::rowvec freq = arma::sum(X, 0) /
                      (2.0 * static_cast<double>(X.n_rows));
  const arma::uvec rare = arma::find(freq < rare_freq);
  const arma::uvec common = arma::find(freq >= rare_freq);

  arma::mat Xnew;
  if (rare.n_elem > 0) {
    arma::mat collapse = arma::sum(X.cols(rare), 1);
    if (common.n_elem > 0) {
      Xnew = arma::join_horiz(X.cols(common), collapse);
    } else {
      Xnew = collapse;
    }
  } else {
    Xnew = X;
  }

  if (Xnew.n_cols == 1) {
    arma::vec collapsed = arma::sum(Xnew, 1);
    collapsed.transform(
        [](double value) { return value > 0.0 ? 1.0 : 0.0; });

    const double total_alt = arma::accu(collapsed);
    const double freq_mutated =
        total_alt / static_cast<double>(collapsed.n_elem);

    const double ncase = arma::accu(Y);
    const double ncont = static_cast<double>(Y.n_elem) - ncase;

    const double case_alt = arma::dot(collapsed, Y);
    const double cont_alt = total_alt - case_alt;
    const double case_ref = ncase - case_alt;
    const double cont_ref = ncont - cont_alt;

    const double case_alt_exp = ncase * freq_mutated;
    const double case_ref_exp = ncase * (1.0 - freq_mutated);
    const double cont_alt_exp = ncont * freq_mutated;
    const double cont_ref_exp = ncont * (1.0 - freq_mutated);

    double stat = std::pow(case_alt - case_alt_exp, 2) / case_alt_exp +
                  std::pow(case_ref - case_ref_exp, 2) / case_ref_exp +
                  std::pow(cont_alt - cont_alt_exp, 2) / cont_alt_exp +
                  std::pow(cont_ref - cont_ref_exp, 2) / cont_ref_exp;

    if (nperm == 0) {
      boost::math::chi_squared chisq(1);
      if (stat < 0.0) {
        stat = std::numeric_limits<double>::epsilon();
      }
      return boost::math::cdf(boost::math::complement(chisq, stat));
    }

    return stat;
  }

  const double n = arma::accu(Xnew);

  arma::rowvec case_counts = Y.t() * Xnew;
  arma::rowvec control_counts = (1.0 - Y).t() * Xnew;

  const double case_total = arma::accu(case_counts);
  const double control_total = arma::accu(control_counts);

  arma::rowvec variant_totals = case_counts + control_counts;

  arma::rowvec case_expected = case_total * variant_totals / n;
  arma::rowvec control_expected = control_total * variant_totals / n;

  arma::rowvec case_chi =
      arma::pow(case_counts - case_expected, 2) / case_expected;
  arma::rowvec control_chi =
      arma::pow(control_counts - control_expected, 2) / control_expected;

  if (nperm == 0) {
    const int df = static_cast<int>(case_counts.n_elem) - 1;
    boost::math::chi_squared chisq(df);
    double stat = arma::accu(case_chi + control_chi);
    if (stat < 0.0) {
      stat = std::numeric_limits<double>::epsilon();
    }
    return boost::math::cdf(boost::math::complement(chisq, stat));
  }

  return arma::accu(case_chi + control_chi);
}
} // namespace

TEST_CASE("Data Construction & Methods") {
  std::stringstream test_data;
  std::stringstream test_cov;
  std::stringstream test_casm;
  std::stringstream test_ped;
  std::string header;

  TaskParams tp;

  tp.maf = 0.5;
  tp.mac = std::numeric_limits<arma::uword>::max();
  tp.method = "VAAST";
  tp.optimizer = "irls";
  tp.no_weights = false;
  tp.aaf_filter = false;
  // File paths and option status
  uint32_t pathbufsize = 1000;
  char pathbuf[pathbufsize];
  for (int i = 0; i < pathbufsize; i++) {
    pathbuf[i] = '\0';
  }
#ifdef __APPLE__
  int ret = _NSGetExecutablePath(pathbuf, &pathbufsize);
#else
  ssize_t len = readlink("/proc/self/exe", pathbuf, 1000);
#endif
  tp.program_path = pathbuf;
  int i = strlen(pathbuf);
  for (; i >= 0; i--) {
    if (pathbuf[i] == '/') {
      break;
    }
  }
  tp.program_directory = std::string(pathbuf).substr(0, i + 1);
  Filter filter(tp.program_directory + "../filter/filter_whitelist.csv");

  // Header
  test_data << "Chr\tStart\tEnd\tRef\tAlt\tType\tGenes\tTranscripts\tRegion\tFu"
               "nction\tAnnotation(c.change:p.change)"
               "\tcase1\tcase2\tcontrol1\tcontrol2\n";
  header = "Chr\tStart\tEnd\tRef\tAlt\tType\tGenes\tTranscripts\tRegion\tFu"
           "nction\tAnnotation(c.change:p.change)"
           "\tcase1\tcase2\tcontrol1\tcontrol2\n";
  // First Transcript
  test_data << "chr1\t1\t1\tA\tC\tSNV\ttest_gene\ttest_"
               "transcript1\tcoding\tnonsynonymous SNV\t.\t0101\n";
  test_data << "chr1\t25\t25\tA\tG\tSNV\ttest_gene\ttest_"
               "transcript1\tcoding\tnonsynonymous SNV\t.\t1111\n";
  test_data << "chr1\t27\t38\tA\tATTACAGATT\tinsertion\ttest_gene\ttest_"
               "transcript1\tcoding\tframeshift\t.\t0222\n";
  test_data << "chr1\t55\t55\tT\tG\tSNV\ttest_gene\ttest_"
               "transcript1\texon\tnonsynonymous SNV\t.\t2001\n";
  // Second Transcript
  test_data << "chr1\t1\t1\tA\tC\tSNV\ttest_gene\ttest_"
               "transcript2\tcoding\tnonsynonymous SNV\t.\t0101\n";
  test_data << "chr1\t25\t25\tA\tG\tSNV\ttest_gene\ttest_"
               "transcript2\tcoding\tnonsynonymous SNV\t.\t1101\n";
  test_data << "chr1\t27\t38\tA\tATTACAGATT\tinsertion\ttest_gene\ttest_"
               "transcript2\tcoding\tframeshift\t.\t0222\n";
  test_data << "chr1\t55\t55\tT\tG\tSNV\ttest_gene\ttest_"
               "transcript2\texon\tnonsynonymous SNV\t.\t2001\n";

  // Phenotypes and Covariates
  test_cov << "control2\t0\t1.5\t1.5\n";
  test_cov << "case2\t1\t0.5\t0.5\n";
  test_cov << "control1\t0\t0.5\t1\n";
  test_cov << "case1\t1\t1\t0.5\n";

  test_ped << "control1\tcontrol1\t0\t0\t0\t1\n";
  test_ped << "control2\tcontrol2\t0\t0\t0\t1\n";
  test_ped << "case1\tcase1\t0\t0\t0\t2\n";
  test_ped << "case2\tcase2\t0\t0\t0\t2\n";

  // Test CASM weights
  // CASM format is chr, start, end, ref, alt, type, gene, transcript, weight
  test_casm << "chr1\t1\t1\tA\tC\tSNV\ttest_gene\ttest_transcript1\t0.5\n";
  test_casm << "chr1\t25\t25\tA\tG\tSNV\ttest_gene\ttest_transcript1\t1.5\n";
  test_casm << "chr1\t27\t38\tA\tATTACAGATT\tinsertion\ttest_gene\ttest_"
               "transcript1\t3.0\n";
  test_casm << "chr1\t55\t55\tT\tG\tSNV\ttest_gene\ttest_transcript1\t0.5\n";
  test_casm << "chr1\t1\t1\tA\tC\tSNV\ttest_gene\ttest_transcript2\t0.5\n";
  test_casm << "chr1\t25\t25\tA\tG\tSNV\ttest_gene\ttest_transcript2\t1.5\n";
  test_casm << "chr1\t27\t38\tA\tATTACAGATT\tinsertion\ttest_gene\ttest_"
               "transcript2\t3.0\n";
  test_casm << "chr1\t55\t55\tT\tG\tSNV\ttest_gene\ttest_transcript2\t0.5\n";

  Covariates cov(test_ped, test_cov, tp);
  cov.sort_covariates(header);
  // Shared pointer for covariates
  std::shared_ptr<Covariates> cov_ptr = std::make_shared<Covariates>(cov);

  // Don't log values
  Weight casm(test_casm);

  // Variant Counts
  std::unordered_map<std::string, arma::uword> nvariants{
      {"test_transcript1", 4}, {"test_transcript2", 4}};

  Gene gene(test_data, cov_ptr, cov.get_nsamples(), nvariants, casm, tp,
            filter);

  SECTION("Data Construction") {

    arma::vec phenotypes{1, 1, 0, 0};
    // Includes a row for the intercept.
    arma::mat covariates{{1.0, 1.0, 1.0, 0.5},
                         {1.0, 1.0, 0.5, 0.5},
                         {1.0, 0.0, 0.5, 1.0},
                         {1.0, 0.0, 1.5, 1.5}};

    arma::mat transcript1{{0.0, 1.0, 2.0, 2.0},
                          {1.0, 1.0, 0.0, 0.0},
                          {0.0, 1.0, 0.0, 0.0}, // Swapped for minor allele
                          {1.0, 1.0, 0.0, 1.0}};
    arma::mat transcript2{{0.0, 1.0, 2.0, 2.0},
                          {1.0, 1.0, 0.0, 0.0}, // Removed missing data
                          {0.0, 0.0, 0.0, 0.0}, // Swapped for minor allele
                          {1.0, 1.0, 0.0, 1.0}};

    arma::vec weights{0.5, 1.5, 3.0, 0.5};

    REQUIRE(cov.get_ncases() == 2);
    REQUIRE(cov.get_nsamples() == 4);
    REQUIRE(arma::all(cov.get_phenotype_vector() == phenotypes));
    REQUIRE(arma::all(arma::all(cov.get_covariate_matrix() == covariates)));

    REQUIRE(gene.gene_name == "test_gene");
    REQUIRE(gene.get_transcripts()[0] == "test_transcript1");
    REQUIRE(gene.get_transcripts()[1] == "test_transcript2");
    REQUIRE(gene.get_transcripts().size() == 2);

    REQUIRE(arma::all(arma::all(arma::mat(gene.genotypes["test_transcript1"]) ==
                                transcript1)));
    REQUIRE(arma::all(arma::all(arma::mat(gene.genotypes["test_transcript2"]) ==
                                transcript2)));
    REQUIRE(arma::all(gene.weights["test_transcript1"] == weights));
  }

  SECTION("BED") {
    std::stringstream test_bed;

    // Test BED file with 5 columns for
    // chromosome, start position, end position, reference allele, alternate
    // allele
    test_bed << "chr1\t1\t1\tA\tT\n";
    test_bed << "chr1\t25\t25\tA\tT\n";
    test_bed << "chr1\t27\t38\tATTACAGATT\tA\n";
    test_bed << "chr1\t55\t55\tA\tT\n";
    test_bed << "chr6\t36867370\t36867370\tA\tT\n";
    test_bed << "chr6\t36867371\t36867371\tA\tT\n";

    Bed bed(test_bed);

    REQUIRE(!bed.empty());
    // Should return the number of variants in the BED file
    REQUIRE(bed.size() == 6);
    // Check for variant as a comma separated string
    REQUIRE(bed.check_variant("chr1,1,1,A,T"));
    REQUIRE(bed.check_variant("chr1,25,25,A,T"));
    REQUIRE(bed.check_variant("chr1,27,38,ATTACAGATT,A"));
    REQUIRE(bed.check_variant("chr1,55,55,A,T"));
    REQUIRE(bed.check_variant("chr6,36867370,36867370,A,T"));
    REQUIRE(bed.check_variant("chr6,36867371,36867371,A,T"));
  }

  SECTION("RANK FUNCTION") {
    arma::vec test_vec{3, 5, 5, 5, 5, 8};
    arma::vec correct{1, 3.5, 3.5, 3.5, 3.5, 6};

    std::cerr << rank(test_vec, "descend", 0);
    REQUIRE(arma::all(rank(test_vec, "ascend", 0) == correct));
    INFO(rank(test_vec, "descend", 0).t());
    REQUIRE(arma::all(rank(test_vec, "descend", 0) == arma::reverse(correct)));

    test_vec = {3, 5, 5, 5, 7, 8};
    correct = {1, 3, 3, 3, 5, 6};
    arma::vec desc_correct = {6, 4, 4, 4, 2, 1};

    REQUIRE(arma::all(rank(test_vec, "ascend", 0) == correct));
    REQUIRE(arma::all(rank(test_vec, "descend", 0) == desc_correct));
  }

  SECTION("CALPHA") {
    tp.method = "CALPHA";
    Methods methods(tp, cov_ptr);

    double calpha =
        methods.CALPHA(gene, cov.get_phenotype_vector(), "test_transcript1");
    REQUIRE(calpha == Approx(-2.000000));

    calpha =
        methods.CALPHA(gene, cov.get_phenotype_vector(), "test_transcript2");
    REQUIRE(calpha == Approx(-1.500000));
  }

  SECTION("SKAT") {
    tp.method = "SKAT";
    tp.kernel = "Linear";
    tp.qtl = false;
    Methods methods(tp, cov_ptr);

    double skat = methods.SKAT(gene, cov.get_phenotype_vector(),
                               "test_transcript1", 1, 25, false, false, false);
    // REQUIRE(skat == Approx(0.625000));
    REQUIRE(skat == Approx(0.0));

    skat = methods.SKAT(gene, cov.get_phenotype_vector(), "test_transcript2", 1,
                        25, false, false, false);
    REQUIRE(skat == Approx(0.0));
  }

  SECTION("VAAST") {
    tp.method = "VAAST";
    tp.group_size = 0;
    Methods methods(tp, cov_ptr);

    double test_vaast =
        methods.VAAST(gene, cov.get_phenotype_vector(), "test_transcript1", 2,
                      0, false, false, 0.5, false);
    REQUIRE(test_vaast == Approx(3.6494094468));

    test_vaast =
        methods.VAAST(gene, cov.get_phenotype_vector(), "test_transcript2", 2,
                      0, false, false, 0.5, false);
    REQUIRE(test_vaast == Approx(3.6494094468));
  }

  SECTION("VT") {
    tp.method = "VT";
    Methods methods(tp, cov_ptr);

    double vt =
        methods.VT(gene, cov.get_phenotype_vector(), "test_transcript1");
    REQUIRE(vt == Approx(0.452267));

    vt = methods.VT(gene, cov.get_phenotype_vector(), "test_transcript2");
    REQUIRE(vt == Approx(0.5345225));
  }

  SECTION("WSS") {
    tp.method = "WSS";
    Methods methods(tp, cov_ptr);

    double wss =
        methods.WSS(gene, cov.get_phenotype_vector(), "test_transcript1");
    REQUIRE(wss == Approx(14.62902));

    wss = methods.WSS(gene, cov.get_phenotype_vector(), "test_transcript2");
    REQUIRE(wss == Approx(14.7115));
  }

  SECTION("CMC sparse partitioning aligns with dense reference") {
    const arma::mat dense = arma::mat(gene.genotypes["test_transcript1"]);
    const arma::vec phen = cov.get_phenotype_vector();

    auto expect_equal = [](double lhs, double rhs) {
      if (std::isnan(rhs)) {
        REQUIRE(std::isnan(lhs));
      } else {
        REQUIRE(lhs == Approx(rhs).epsilon(1e-10));
      }
    };

    TaskParams cmc_multi{};
    cmc_multi.method = "CMC";
    cmc_multi.nperm = 0;
    Methods multi_methods(cmc_multi, cov_ptr);
    const double multi_threshold = 0.4;
    arma::vec phen_multi = phen;
    const double sparse_multi =
        multi_methods.CMC(gene, phen_multi, "test_transcript1", multi_threshold);
    const double dense_multi =
        cmc_dense_reference(dense, phen, multi_threshold, false, 0);
    expect_equal(sparse_multi, dense_multi);

    TaskParams cmc_switch{};
    cmc_switch.method = "CMC";
    cmc_switch.nperm = 5;
    Methods switch_methods(cmc_switch, cov_ptr);
    const double switch_threshold = 0.6;
    arma::vec phen_switch = phen;
    const double sparse_switch = switch_methods.CMC(
        gene, phen_switch, "test_transcript1", switch_threshold);
    const double dense_switch =
        cmc_dense_reference(dense, phen, switch_threshold, false, 5);
    expect_equal(sparse_switch, dense_switch);

    TaskParams cmc_hot{};
    cmc_hot.method = "CMC";
    cmc_hot.hotellings = true;
    cmc_hot.nperm = 0;
    Methods hot_methods(cmc_hot, cov_ptr);
    const double hot_threshold = 0.4;
    arma::vec phen_hot = phen;
    const double sparse_hot =
        hot_methods.CMC(gene, phen_hot, "test_transcript1", hot_threshold);
    const double dense_hot =
        cmc_dense_reference(dense, phen, hot_threshold, true, 0);
    expect_equal(sparse_hot, dense_hot);

    cmc_hot.nperm = 11;
    Methods hot_stat_methods(cmc_hot, cov_ptr);
    arma::vec phen_hot_stat = phen;
    const double sparse_hot_stat = hot_stat_methods.CMC(
        gene, phen_hot_stat, "test_transcript1", hot_threshold);
    const double dense_hot_stat =
        cmc_dense_reference(dense, phen, hot_threshold, true, 11);
    expect_equal(sparse_hot_stat, dense_hot_stat);
  }

  SECTION("CMC") {
    arma::sp_mat sim_data(arma::mat{
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 1, 0, 0, 1, 0, 0, 0, 0, 0},
        {0, 1, 0, 0, 0, 0, 0, 0, 0, 1}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 1, 0, 0, 0, 0, 0, 0, 0, 0}, {2, 0, 0, 0, 0, 0, 0, 1, 0, 0},
        {2, 0, 1, 0, 1, 0, 0, 0, 0, 0}, {1, 1, 0, 1, 1, 0, 0, 0, 0, 1},
        {0, 0, 1, 1, 0, 0, 0, 1, 0, 1}, {0, 1, 0, 0, 1, 0, 0, 0, 0, 1},
        {0, 1, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 1, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 1, 0, 0, 2, 0, 0, 0, 0, 0}, {0, 1, 0, 0, 0, 0, 0, 0, 0, 1},
        {0, 0, 0, 0, 0, 0, 1, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 1}, {0, 0, 1, 0, 1, 1, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 1, 0, 0, 0, 0, 1, 0, 0, 0}, {0, 0, 1, 0, 1, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 1, 0, 0, 0, 0, 0}, {1, 0, 1, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 1, 1, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 1, 0, 0, 0, 0, 0, 0, 0, 0}, {1, 0, 1, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 1, 0, 0, 0, 0, 1}, {1, 0, 0, 0, 1, 0, 0, 0, 0, 0},
        {1, 1, 1, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 1, 1, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 1, 1, 0, 0, 1}, {1, 0, 0, 0, 1, 1, 0, 0, 0, 0},
        {0, 0, 0, 1, 1, 0, 1, 0, 0, 0}, {1, 0, 0, 1, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 1, 0, 1, 0, 0, 0},
        {0, 0, 0, 0, 1, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 2},
        {0, 0, 0, 0, 1, 0, 0, 1, 0, 0}, {1, 0, 1, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 1, 0, 0, 0, 0, 0}, {2, 0, 1, 1, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 1, 1, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 1, 0, 1, 0, 0, 0, 0, 0, 0}, {0, 1, 0, 0, 1, 0, 0, 0, 0, 0},
        {0, 0, 0, 1, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 1}, {0, 1, 0, 0, 0, 1, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 1, 0, 0, 0}, {1, 0, 1, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 1}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 1, 0, 1, 0, 0, 0}, {0, 2, 1, 0, 0, 0, 0, 0, 0, 0},
        {1, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 1, 0, 1, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 1, 0, 1},
        {0, 1, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 1, 1, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 1, 0, 1, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 1, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 1}, {0, 0, 1, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 1, 0, 1, 0, 0, 0, 1}, {1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 1, 0, 0, 0, 0, 1, 0, 0, 1},
        {0, 0, 0, 1, 0, 0, 0, 0, 0, 1}, {0, 0, 0, 1, 0, 1, 2, 0, 0, 1},
        {0, 1, 0, 1, 0, 0, 0, 0, 0, 0}, {1, 0, 0, 0, 0, 1, 0, 0, 0, 0},
        {0, 0, 1, 0, 0, 0, 0, 0, 0, 1}, {1, 0, 0, 0, 0, 0, 1, 0, 0, 0},
        {0, 1, 0, 0, 1, 0, 0, 1, 0, 0}, {0, 0, 2, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 1, 0, 0, 1, 0, 0, 0}, {0, 1, 0, 0, 1, 0, 0, 0, 0, 0},
        {1, 0, 0, 0, 1, 0, 0, 0, 0, 0}, {0, 1, 0, 0, 1, 0, 0, 0, 0, 0},
        {1, 1, 1, 0, 1, 0, 0, 0, 0, 1}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 1, 1, 0, 0, 1, 0, 0}, {0, 1, 0, 1, 1, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 1, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 1}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 1, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 1, 0, 0, 1, 0, 1, 1, 0, 0}, {0, 0, 0, 0, 1, 0, 0, 0, 0, 1},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {1, 0, 0, 0, 0, 0, 0, 1, 0, 2},
        {0, 0, 1, 0, 0, 0, 0, 0, 0, 2}, {0, 0, 0, 0, 1, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 1, 1, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 1, 0, 1, 0, 1, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 1, 0, 0, 0, 0, 0, 0, 0}, {2, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {1, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 1, 0, 0, 1, 0, 0, 0, 0, 0},
        {0, 1, 0, 0, 0, 0, 0, 0, 0, 1}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 1},
        {0, 0, 0, 1, 1, 1, 0, 0, 0, 0}, {0, 0, 0, 0, 1, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 1, 0, 0, 0, 0, 0}, {1, 0, 0, 0, 0, 0, 0, 0, 0, 1},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 1, 1, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 2}, {1, 0, 0, 0, 0, 0, 1, 0, 0, 0},
        {2, 0, 0, 0, 0, 0, 0, 0, 0, 1}, {0, 0, 0, 0, 1, 0, 0, 0, 0, 1},
        {0, 0, 0, 0, 1, 1, 0, 0, 0, 0}, {0, 0, 0, 1, 0, 0, 0, 0, 0, 1},
        {1, 0, 1, 1, 0, 1, 0, 0, 0, 1}, {0, 0, 0, 1, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 1, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {1, 1, 1, 1, 1, 0, 0, 0, 0, 0}, {0, 0, 1, 0, 0, 0, 0, 0, 0, 1},
        {0, 0, 0, 0, 2, 0, 0, 0, 0, 1}, {0, 1, 0, 1, 1, 0, 0, 0, 0, 0},
        {0, 0, 0, 1, 0, 0, 0, 0, 0, 0}, {1, 1, 1, 0, 1, 0, 0, 0, 0, 1},
        {1, 1, 0, 0, 1, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 1, 2, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 1, 0, 0, 0, 0},
        {0, 1, 0, 0, 1, 0, 0, 0, 0, 0}, {1, 0, 0, 1, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 2}, {1, 1, 0, 0, 1, 0, 0, 0, 0, 1},
        {0, 0, 0, 0, 1, 0, 0, 0, 0, 1}, {0, 0, 0, 0, 0, 1, 0, 0, 0, 0},
        {0, 0, 0, 0, 1, 0, 1, 0, 0, 0}, {0, 0, 1, 0, 0, 0, 0, 0, 0, 1},
        {0, 1, 0, 0, 0, 0, 0, 0, 0, 1}, {0, 0, 1, 0, 1, 0, 0, 0, 0, 1},
        {0, 0, 1, 1, 0, 0, 0, 0, 0, 0}, {1, 1, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 1, 0, 0, 0, 0, 0, 0, 1}, {1, 0, 1, 1, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 1, 0, 1, 2, 0, 0, 1}, {0, 1, 0, 0, 0, 0, 0, 0, 0, 1},
        {0, 0, 1, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 1, 0, 0, 0, 0},
        {0, 0, 1, 0, 0, 0, 0, 0, 0, 0}, {0, 1, 0, 1, 0, 0, 0, 0, 0, 1},
        {1, 0, 0, 1, 0, 0, 0, 0, 0, 0}, {0, 1, 0, 1, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 1, 0, 1},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 1, 0, 0, 1, 0, 1, 0, 0, 0},
        {0, 0, 0, 1, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 1},
        {1, 1, 1, 1, 0, 0, 0, 0, 0, 1}, {0, 0, 0, 0, 0, 1, 0, 0, 0, 0},
        {1, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {1, 0, 0, 1, 0, 0, 0, 0, 0, 1}, {0, 0, 0, 0, 1, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 1, 0, 0, 0}, {0, 0, 1, 0, 0, 0, 1, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 1, 0, 0, 0}, {0, 0, 1, 0, 0, 0, 0, 0, 0, 1},
        {0, 0, 1, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 2, 0, 0, 0, 0, 0}, {0, 1, 0, 0, 1, 0, 0, 0, 0, 1},
        {0, 1, 0, 0, 0, 0, 1, 0, 0, 1}, {0, 0, 0, 1, 1, 0, 0, 0, 0, 0},
        {1, 1, 1, 0, 1, 1, 0, 0, 0, 0}, {1, 0, 0, 0, 0, 1, 0, 0, 0, 0},
        {0, 0, 0, 0, 1, 0, 0, 0, 0, 0}, {0, 1, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 1, 0, 0, 0},
        {1, 0, 2, 0, 0, 0, 0, 0, 0, 1}, {0, 0, 0, 0, 1, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 1, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 1, 0, 0, 0, 0, 1},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 1, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {1, 0, 0, 0, 0, 0, 0, 0, 0, 1},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 1, 0, 0, 0, 0, 0, 0, 0, 0},
        {1, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 2, 0, 1, 0, 0, 0, 0, 0},
        {1, 0, 0, 0, 0, 0, 0, 0, 0, 1}, {0, 0, 0, 1, 0, 0, 0, 0, 0, 0},
        {0, 1, 1, 0, 0, 1, 0, 0, 0, 1}, {0, 0, 0, 0, 0, 1, 0, 0, 0, 0},
        {0, 0, 0, 1, 0, 0, 0, 0, 0, 0}, {0, 1, 0, 0, 0, 0, 0, 0, 0, 0},
        {1, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 1, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 1}, {0, 0, 1, 0, 1, 0, 0, 0, 0, 1},
        {1, 1, 1, 0, 0, 0, 0, 0, 0, 0}, {1, 1, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 1, 0, 0},
        {1, 0, 0, 0, 1, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 1, 0, 0, 0},
        {0, 1, 0, 0, 1, 0, 0, 0, 0, 0}, {0, 0, 1, 0, 1, 0, 0, 0, 0, 0},
        {1, 0, 1, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 2, 0, 0, 0, 0, 1},
        {1, 1, 1, 1, 1, 1, 1, 0, 0, 0}, {0, 0, 0, 0, 1, 0, 0, 0, 0, 1},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 2}, {0, 0, 0, 0, 1, 0, 0, 0, 0, 1},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {1, 0, 0, 0, 0, 1, 0, 0, 0, 0},
        {0, 0, 0, 0, 1, 0, 1, 0, 0, 0}, {0, 0, 0, 1, 0, 0, 0, 0, 0, 0},
        {0, 0, 1, 0, 2, 0, 0, 0, 0, 0}, {0, 0, 0, 1, 0, 0, 0, 0, 0, 2},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 1, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 1, 1, 0, 1},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {1, 0, 0, 0, 0, 1, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 1, 0, 0, 0}, {1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 1, 0, 0, 0, 1, 0, 0, 0, 1}, {0, 0, 0, 0, 0, 0, 0, 1, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 1}, {1, 0, 0, 0, 0, 0, 1, 0, 0, 0},
        {0, 0, 1, 0, 0, 0, 0, 0, 0, 1}, {0, 1, 0, 1, 0, 1, 0, 0, 0, 0},
        {0, 0, 0, 1, 1, 0, 1, 0, 0, 0}, {0, 1, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 1, 0, 0, 0, 0, 0, 0, 2},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 1, 0, 0, 0},
        {0, 0, 0, 1, 0, 0, 0, 0, 0, 1}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 1, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 1, 0, 0, 1, 0, 0, 0, 0},
        {1, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 1, 0, 1, 0, 0, 1, 0, 0},
        {1, 0, 0, 0, 2, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 1, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 1, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 1, 0, 0, 1, 0, 0, 1, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 1},
        {0, 0, 0, 1, 0, 0, 1, 0, 0, 0}, {0, 0, 0, 0, 1, 1, 0, 0, 0, 0},
        {0, 1, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 1, 0, 2, 0, 0, 0},
        {1, 0, 0, 1, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 1, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 1, 0, 0, 0, 0, 0, 0, 0, 2},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 1},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 1}, {0, 0, 0, 0, 1, 0, 0, 0, 0, 0},
        {0, 1, 0, 0, 2, 0, 0, 0, 0, 1}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 1, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 1, 0, 0, 0, 0, 1},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {1, 1, 0, 0, 1, 0, 0, 0, 0, 0},
        {0, 1, 0, 0, 0, 2, 0, 0, 0, 1}, {0, 0, 0, 0, 0, 1, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 1},
        {0, 0, 0, 1, 1, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 1, 0, 0, 0, 0, 1},
        {0, 0, 1, 0, 0, 0, 2, 0, 0, 1}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {1, 1, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {1, 0, 1, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 1},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 1},
        {0, 0, 0, 1, 1, 0, 1, 0, 0, 0}, {1, 2, 1, 1, 1, 0, 0, 0, 0, 1},
        {0, 0, 1, 1, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 1, 0, 0, 0, 0, 0},
        {0, 0, 0, 1, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 1, 0, 0, 0, 0, 0, 0},
        {0, 0, 1, 0, 1, 0, 1, 0, 0, 0}, {1, 1, 0, 1, 0, 1, 0, 0, 0, 1},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 1, 0, 0, 1, 0, 0, 0},
        {0, 0, 0, 0, 1, 0, 0, 0, 0, 0}, {0, 0, 1, 0, 0, 0, 1, 0, 0, 0},
        {0, 0, 1, 2, 1, 0, 0, 0, 0, 1}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 1},
        {1, 0, 0, 1, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 1, 0, 1, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 1, 0, 0, 0, 0, 0, 1},
        {1, 0, 0, 0, 0, 0, 1, 0, 0, 1}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 1, 0, 1, 1, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 1, 0, 0, 1}, {0, 0, 0, 1, 0, 0, 0, 1, 0, 0},
        {0, 0, 1, 0, 1, 0, 0, 0, 0, 0}, {0, 1, 0, 0, 0, 0, 0, 0, 0, 0},
        {1, 0, 0, 0, 0, 0, 1, 0, 0, 0}, {0, 1, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 1, 0, 0, 0, 0, 1, 0, 0, 0}, {0, 0, 1, 0, 1, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 1, 0, 0, 0, 0, 1}, {0, 0, 1, 1, 2, 0, 0, 0, 0, 1},
        {0, 0, 0, 0, 1, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 1, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 2}, {0, 0, 0, 1, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 1, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 2, 0, 0, 1, 0, 0, 0, 0},
        {0, 0, 0, 1, 0, 1, 0, 0, 0, 2}, {0, 0, 0, 1, 1, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 1, 0, 1, 0, 0, 0, 0},
        {0, 0, 1, 0, 0, 0, 0, 0, 0, 1}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {1, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 1},
        {1, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 1},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 1, 0, 0, 0, 0, 0}, {0, 1, 0, 1, 0, 0, 1, 0, 0, 0},
        {0, 0, 0, 0, 1, 0, 1, 0, 0, 0}, {1, 1, 0, 0, 0, 1, 0, 0, 0, 2},
        {0, 0, 0, 0, 0, 1, 0, 0, 0, 1}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 1, 0, 0, 0, 0}, {0, 0, 1, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 2, 0, 2, 0, 0, 0}, {0, 0, 1, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 1, 0, 0, 0, 0, 1},
        {1, 0, 0, 0, 0, 0, 0, 1, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 1, 0, 0, 0, 0, 0, 0, 0, 2}, {2, 0, 0, 1, 0, 1, 0, 0, 0, 0},
        {0, 1, 0, 0, 0, 0, 0, 0, 0, 0}, {1, 0, 1, 0, 0, 0, 0, 0, 0, 0},
        {0, 1, 0, 0, 0, 0, 1, 0, 0, 0}, {1, 0, 0, 0, 1, 0, 0, 0, 0, 1},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 1, 1, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 1, 1, 0, 0, 0, 0, 0}, {0, 1, 1, 0, 0, 0, 0, 0, 0, 0},
        {1, 0, 0, 0, 0, 0, 0, 1, 0, 0}, {2, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 1}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 1, 1, 0, 0}, {0, 1, 0, 0, 1, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 1, 0, 0, 1, 0, 0}, {1, 0, 1, 0, 0, 0, 0, 0, 0, 1},
        {0, 0, 0, 2, 0, 0, 0, 0, 0, 0}, {1, 1, 0, 1, 0, 0, 0, 0, 0, 1},
        {0, 0, 0, 1, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {1, 0, 0, 0, 1, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 1},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 1}, {1, 0, 1, 1, 1, 0, 0, 0, 0, 0},
        {0, 0, 0, 1, 0, 0, 0, 0, 0, 0}, {1, 0, 0, 0, 0, 1, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {1, 1, 0, 1, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 1, 0, 0, 0, 0, 1}, {1, 0, 0, 0, 0, 0, 1, 0, 0, 0},
        {0, 0, 2, 1, 1, 1, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {1, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 1, 1, 0, 0, 0, 0, 1},
        {2, 1, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 1, 0, 0, 1, 0, 1, 0, 0, 1},
        {0, 1, 1, 1, 0, 0, 0, 0, 0, 0}, {0, 1, 0, 1, 0, 0, 0, 0, 0, 0},
        {1, 0, 1, 0, 0, 0, 0, 0, 0, 0}, {0, 1, 1, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 1, 0, 0, 0, 1, 0, 0, 0}, {0, 1, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 2, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 1},
        {0, 0, 0, 0, 1, 0, 0, 1, 0, 0}, {0, 1, 0, 1, 0, 0, 0, 0, 0, 0},
        {0, 2, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 1, 0, 1, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 1, 0, 0, 0, 0, 1, 0, 1},
        {0, 0, 1, 0, 0, 0, 0, 0, 0, 1}, {0, 0, 1, 0, 0, 0, 0, 0, 0, 1},
        {0, 0, 1, 0, 1, 0, 0, 0, 0, 0}, {0, 1, 0, 0, 0, 1, 0, 0, 0, 0},
        {0, 1, 0, 0, 1, 0, 0, 1, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 1, 0, 0, 0, 0, 1, 1, 0, 0}, {1, 1, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 1, 0, 1, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 1, 2, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 1, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 1},
        {0, 1, 0, 0, 1, 1, 0, 0, 0, 0}, {0, 0, 1, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 1, 0, 0}, {1, 0, 0, 0, 1, 1, 0, 0, 0, 0},
        {0, 0, 0, 0, 1, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 1, 0, 1, 0, 0, 0},
        {1, 0, 0, 0, 0, 0, 0, 0, 0, 2}, {0, 0, 1, 0, 0, 0, 1, 2, 0, 0},
        {0, 0, 0, 0, 0, 0, 1, 1, 0, 0}, {0, 1, 1, 0, 1, 0, 0, 0, 0, 0},
        {0, 0, 0, 1, 0, 0, 0, 0, 0, 1}, {0, 0, 1, 1, 0, 0, 0, 0, 0, 0},
        {0, 1, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 2, 1, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 2}, {0, 0, 1, 0, 0, 0, 0, 0, 0, 0},
        {0, 1, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 1, 0, 0, 0, 0, 1},
        {1, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 1},
        {0, 1, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 1, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 2, 0, 1, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 1, 0, 0, 0, 0, 1}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 2},
        {0, 0, 1, 0, 1, 0, 0, 0, 0, 2}, {0, 2, 1, 1, 0, 0, 0, 0, 0, 1},
        {0, 0, 1, 0, 0, 0, 0, 0, 0, 1}, {0, 0, 0, 1, 1, 1, 0, 0, 0, 0},
        {0, 0, 1, 0, 0, 0, 0, 0, 0, 1}, {0, 0, 1, 1, 1, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 2, 0, 1, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 1}, {0, 1, 0, 1, 1, 0, 0, 0, 0, 0},
        {1, 1, 0, 0, 0, 0, 0, 1, 0, 1}, {0, 0, 1, 1, 1, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 1}, {0, 0, 1, 0, 0, 0, 1, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 1, 0, 1}, {0, 0, 0, 0, 1, 0, 0, 0, 0, 1},
        {0, 0, 0, 1, 1, 0, 0, 0, 0, 1}, {0, 1, 1, 1, 1, 0, 0, 0, 0, 0},
        {0, 2, 1, 1, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 1, 1, 1, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 2, 0, 0, 0, 0, 1},
        {0, 0, 0, 0, 1, 0, 1, 0, 0, 0}, {0, 0, 0, 0, 1, 0, 0, 1, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 1, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 1, 0, 0}, {1, 0, 0, 0, 1, 0, 0, 0, 0, 0},
        {0, 1, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 1, 0, 0, 1, 0, 0, 1},
        {1, 1, 0, 0, 0, 0, 0, 0, 0, 1}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 1, 0, 1, 0, 0, 0, 0, 0},
        {1, 0, 0, 0, 1, 0, 1, 0, 0, 0}, {0, 0, 2, 0, 0, 0, 0, 0, 0, 0},
        {1, 0, 0, 0, 0, 0, 1, 0, 0, 0}, {0, 1, 0, 0, 1, 0, 0, 0, 0, 1},
        {0, 0, 0, 1, 0, 0, 0, 0, 0, 0}, {1, 0, 0, 0, 0, 0, 0, 0, 0, 1},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 2}, {0, 0, 0, 1, 0, 0, 0, 1, 0, 0},
        {0, 0, 1, 0, 0, 0, 0, 1, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 1, 0, 0, 2, 0, 0, 0, 0, 0}, {0, 1, 0, 1, 1, 0, 0, 0, 0, 0},
        {1, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 1, 0, 0, 0, 0, 0}, {1, 1, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 1, 0, 0, 0, 0, 1}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 1},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {1, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 1, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 1, 0, 1, 0, 0, 0, 0, 0},
        {0, 0, 0, 1, 0, 0, 0, 0, 0, 1}, {0, 1, 0, 1, 1, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {1, 0, 0, 0, 1, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 1, 0, 1, 0, 0, 1},
        {0, 0, 0, 1, 0, 0, 1, 0, 0, 0}, {0, 0, 0, 1, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 1, 0, 0, 0, 0, 0, 2}, {0, 1, 0, 1, 1, 0, 1, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 1, 0, 0, 1, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 2, 0, 0, 0, 0, 2},
        {0, 0, 0, 1, 1, 0, 0, 1, 0, 1}, {0, 1, 0, 0, 0, 1, 0, 0, 0, 1},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 1, 1, 0, 1, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 1, 0, 0, 0, 0, 0}, {1, 0, 1, 0, 0, 1, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 1, 0, 0, 0, 0, 0}, {0, 1, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 1, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 1, 0, 0, 0, 0, 0, 0, 2}, {0, 0, 2, 1, 0, 0, 0, 0, 0, 0},
        {0, 1, 1, 0, 0, 0, 1, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 1, 0, 0, 0, 0, 1, 0, 0, 0}, {1, 1, 0, 0, 1, 0, 1, 0, 0, 0},
        {1, 0, 1, 1, 1, 0, 0, 0, 0, 0}, {0, 0, 0, 1, 1, 0, 0, 0, 0, 1},
        {0, 0, 1, 0, 2, 1, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 1, 0, 1, 0, 0, 0, 0, 0}, {0, 1, 0, 1, 2, 0, 0, 0, 0, 0},
        {0, 0, 1, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 1, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 1, 0, 0, 1, 0, 0, 0, 0, 0},
        {0, 0, 1, 0, 0, 0, 0, 1, 0, 0}, {0, 0, 0, 1, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 1, 0, 0, 0, 0, 0, 1}, {0, 0, 2, 1, 0, 0, 0, 0, 0, 1},
        {0, 0, 0, 0, 1, 0, 0, 1, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 1},
        {0, 0, 0, 0, 1, 0, 0, 0, 0, 0}, {1, 1, 0, 0, 0, 0, 0, 0, 0, 1},
        {0, 0, 1, 1, 0, 0, 0, 0, 0, 0}, {1, 0, 1, 1, 1, 0, 1, 0, 0, 0},
        {0, 0, 1, 0, 0, 0, 0, 0, 0, 0}, {0, 1, 1, 0, 1, 0, 1, 0, 0, 1},
        {0, 0, 1, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 1, 0, 0, 1, 0, 0, 0, 0},
        {1, 0, 0, 0, 1, 0, 0, 0, 0, 1}, {0, 1, 0, 0, 0, 0, 0, 0, 0, 0},
        {1, 0, 0, 0, 0, 1, 0, 0, 0, 0}, {0, 1, 0, 1, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 1, 0, 0, 0, 0, 0, 1}, {0, 1, 0, 0, 1, 1, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 1, 0, 1, 0, 0, 0, 0, 0, 1},
        {0, 0, 0, 0, 0, 0, 1, 0, 0, 1}, {0, 0, 1, 0, 1, 0, 1, 1, 0, 0},
        {0, 1, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 1, 0, 0, 0, 0, 0, 0},
        {1, 1, 0, 0, 1, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 1, 0, 0, 0, 0},
        {0, 1, 1, 0, 1, 1, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 1, 0, 1, 1, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 1, 0, 1, 0, 0, 0, 0, 0, 1}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 1, 0, 0, 0, 1, 0, 0, 0}, {0, 0, 0, 1, 0, 0, 0, 0, 0, 0},
        {1, 0, 0, 1, 0, 0, 1, 0, 0, 0}, {0, 0, 1, 1, 0, 0, 0, 0, 0, 0},
        {1, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 1, 0, 0, 1, 0, 0, 0, 0, 0},
        {0, 1, 0, 0, 0, 0, 1, 0, 0, 1}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 1, 0, 0}, {0, 0, 0, 1, 0, 0, 0, 0, 0, 1},
        {0, 1, 0, 0, 0, 0, 1, 0, 0, 0}, {0, 2, 0, 1, 1, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 1, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 1},
        {0, 0, 0, 0, 1, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 1},
        {0, 0, 0, 1, 0, 0, 1, 0, 0, 0}, {0, 0, 0, 0, 1, 0, 0, 0, 0, 1},
        {0, 1, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 1, 0, 0, 1, 0, 0},
        {0, 0, 0, 0, 0, 1, 0, 0, 0, 1}, {0, 0, 1, 0, 0, 0, 0, 1, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 1}, {0, 0, 0, 0, 0, 1, 0, 0, 0, 0},
        {0, 0, 1, 2, 0, 0, 0, 1, 0, 0}, {1, 0, 0, 0, 1, 0, 0, 0, 0, 0},
        {0, 1, 0, 0, 1, 0, 1, 0, 0, 0}, {1, 0, 1, 0, 1, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 2, 1, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 1},
        {0, 0, 0, 0, 0, 0, 0, 1, 0, 0}, {1, 0, 0, 0, 1, 1, 0, 1, 0, 0},
        {1, 0, 0, 0, 0, 0, 0, 0, 0, 2}, {0, 1, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 1, 1, 1, 0, 0, 0, 1, 0, 0}, {0, 0, 1, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 1, 0, 0, 0, 0, 0, 0},
        {0, 0, 1, 0, 0, 0, 0, 0, 0, 0}, {0, 1, 0, 1, 1, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 1}, {0, 0, 0, 0, 0, 1, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {1, 1, 1, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 1, 1},
        {1, 1, 0, 0, 0, 0, 1, 0, 0, 0}, {0, 1, 0, 0, 0, 0, 0, 1, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {1, 0, 1, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 1, 0, 0, 0, 1, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 1},
        {0, 0, 0, 1, 0, 0, 1, 0, 0, 0}, {1, 0, 1, 0, 1, 0, 1, 0, 0, 1},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {1, 0, 0, 0, 0, 1, 0, 0, 0, 1},
        {0, 0, 0, 0, 0, 0, 1, 0, 0, 0}, {0, 1, 0, 0, 1, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 1, 0, 0, 0, 0, 0}, {1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 1, 0, 0, 1}, {0, 0, 0, 0, 1, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 1, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 1, 0, 0, 1},
        {0, 0, 1, 1, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 1, 0, 0, 0, 0, 0, 1},
        {0, 0, 0, 2, 1, 0, 0, 0, 0, 0}, {0, 1, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 1, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 1, 0, 0, 0, 0, 0, 0, 0, 2}, {0, 0, 0, 1, 0, 0, 1, 0, 1, 0},
        {1, 2, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 1, 0, 1, 0, 0, 0, 0, 1},
        {0, 0, 0, 0, 0, 1, 1, 0, 0, 0}, {0, 0, 1, 0, 0, 0, 0, 1, 0, 0},
        {0, 0, 2, 1, 0, 0, 0, 0, 0, 0}, {0, 0, 1, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 2, 0, 0, 0, 0, 0}, {0, 0, 0, 1, 1, 0, 0, 0, 0, 0},
        {1, 1, 0, 0, 0, 0, 0, 1, 0, 0}, {1, 1, 0, 0, 1, 1, 0, 0, 0, 0},
        {0, 0, 0, 0, 1, 0, 0, 0, 0, 1}, {0, 0, 2, 1, 1, 0, 0, 1, 0, 1},
        {0, 0, 0, 1, 0, 1, 0, 0, 0, 0}, {0, 0, 1, 1, 0, 0, 0, 0, 0, 1},
        {0, 0, 0, 0, 0, 0, 1, 0, 0, 1}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 1, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 1, 0, 0, 1, 0, 0, 1, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 1}, {0, 0, 0, 1, 0, 0, 0, 0, 0, 1},
        {0, 0, 1, 1, 0, 0, 0, 1, 0, 1}, {0, 0, 1, 0, 0, 1, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 1, 0}, {0, 0, 0, 1, 0, 0, 0, 0, 0, 0},
        {0, 1, 0, 0, 0, 0, 0, 1, 0, 0}, {0, 0, 1, 0, 0, 0, 0, 0, 0, 0},
        {0, 1, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 1, 1, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 2}, {0, 0, 0, 0, 1, 1, 0, 0, 0, 0},
        {0, 0, 1, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 1, 0, 0, 0, 0, 0, 0, 0},
        {0, 1, 0, 1, 1, 0, 1, 0, 0, 0}, {0, 1, 1, 0, 0, 1, 0, 0, 0, 0},
        {1, 0, 0, 0, 0, 1, 0, 0, 0, 1}, {0, 1, 1, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 2, 0, 2, 0, 0, 1, 0, 0}, {1, 0, 0, 1, 1, 0, 1, 0, 0, 0},
        {1, 0, 1, 0, 0, 0, 0, 0, 0, 1}, {0, 0, 0, 0, 1, 0, 0, 0, 0, 1},
        {0, 0, 0, 0, 1, 0, 1, 0, 0, 0}, {1, 0, 0, 1, 0, 0, 0, 0, 0, 1},
        {0, 0, 1, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 1, 0, 0, 0, 0, 0, 0, 0},
        {0, 1, 0, 0, 1, 0, 0, 0, 0, 1}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 1, 0, 1, 0, 0, 1, 0, 1}, {1, 0, 1, 0, 1, 1, 0, 0, 0, 0},
        {0, 0, 0, 0, 1, 0, 0, 0, 0, 0}, {0, 1, 1, 0, 0, 0, 0, 0, 0, 1},
        {0, 0, 0, 1, 1, 0, 0, 0, 0, 0}, {0, 0, 1, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 1}, {0, 0, 0, 1, 0, 0, 2, 0, 0, 0},
        {0, 0, 0, 0, 1, 0, 0, 0, 1, 0}, {2, 0, 0, 0, 0, 0, 0, 0, 0, 1},
        {0, 1, 0, 0, 0, 0, 0, 0, 0, 1}, {0, 1, 1, 1, 0, 0, 1, 0, 0, 0},
        {0, 1, 0, 0, 0, 0, 0, 0, 0, 0}, {1, 0, 0, 0, 0, 1, 0, 1, 0, 0},
        {0, 0, 1, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 1, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 1, 0}, {0, 1, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 1, 1, 0, 0, 0, 0, 0, 0, 1}, {0, 1, 1, 2, 0, 0, 0, 0, 0, 0},
        {1, 0, 0, 1, 0, 0, 0, 0, 1, 0}, {0, 1, 0, 0, 0, 1, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 1, 0, 1, 0, 0, 0, 1},
        {1, 0, 0, 0, 1, 0, 0, 0, 0, 0}, {1, 0, 0, 0, 1, 0, 0, 1, 0, 1},
        {0, 0, 0, 0, 1, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 1, 0, 1, 0, 0, 0},
        {1, 0, 1, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 1, 0, 0, 0, 0, 0},
        {1, 0, 1, 0, 0, 0, 0, 0, 0, 0}, {1, 1, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 1, 0, 1, 0, 0, 0, 1}, {0, 0, 0, 1, 1, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 1}, {1, 0, 0, 0, 1, 0, 0, 0, 0, 0},
        {0, 0, 0, 1, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 1, 1, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 1, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 1, 0, 0, 0, 0, 0, 0, 1}, {0, 0, 0, 1, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 1, 0, 0, 0, 0, 0, 0, 0, 1},
        {0, 0, 1, 0, 1, 0, 0, 0, 0, 0}, {0, 1, 0, 0, 2, 0, 0, 0, 0, 1},
        {0, 0, 0, 0, 1, 0, 0, 0, 0, 0}, {0, 1, 0, 1, 1, 0, 0, 0, 0, 0},
        {0, 1, 0, 1, 0, 0, 0, 0, 0, 0}, {0, 0, 1, 0, 0, 0, 0, 0, 1, 0},
        {1, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 1, 0, 0, 0, 0, 1},
        {0, 0, 0, 0, 0, 0, 0, 1, 0, 0}, {0, 1, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 1, 0, 1, 0, 0, 0}, {1, 1, 2, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 1, 0, 0, 0, 0, 0, 1}, {1, 0, 1, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 1, 0, 0, 1, 0, 0, 0}, {0, 0, 1, 0, 0, 0, 0, 0, 0, 0},
        {1, 1, 0, 0, 1, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 1},
        {0, 0, 0, 0, 1, 0, 0, 0, 0, 0}, {0, 0, 1, 0, 0, 0, 1, 0, 0, 0},
        {1, 0, 0, 0, 1, 0, 0, 0, 0, 0}, {0, 0, 1, 0, 1, 0, 1, 0, 0, 0},
        {0, 0, 1, 1, 0, 0, 0, 0, 0, 1}, {1, 1, 0, 0, 0, 0, 1, 0, 0, 1},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 1, 0, 0, 1},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 1}, {0, 0, 0, 0, 1, 0, 1, 0, 0, 0},
        {0, 0, 0, 0, 1, 0, 1, 1, 0, 1}, {1, 1, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 1, 0, 0, 1, 1, 0, 0, 0, 0}, {0, 0, 0, 1, 0, 0, 1, 0, 0, 0},
        {0, 1, 0, 2, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {1, 0, 1, 1, 0, 0, 0, 0, 0, 1}, {0, 1, 0, 0, 0, 0, 0, 0, 0, 0},
        {1, 0, 0, 0, 2, 0, 0, 0, 0, 1}, {0, 0, 0, 0, 1, 0, 0, 0, 0, 0},
        {0, 0, 2, 1, 0, 0, 0, 0, 0, 0}, {0, 2, 0, 0, 0, 0, 1, 0, 0, 0},
        {1, 0, 0, 0, 0, 0, 1, 0, 0, 0}, {1, 0, 0, 0, 0, 0, 0, 1, 0, 1},
        {0, 0, 0, 0, 1, 0, 0, 0, 0, 0}, {0, 0, 1, 1, 1, 0, 0, 0, 0, 1},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 1}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 1, 1, 0, 0, 0, 0, 1},
        {1, 0, 0, 0, 1, 1, 0, 0, 0, 0}, {0, 1, 0, 0, 1, 0, 0, 0, 0, 0},
        {1, 1, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 1, 0, 0, 0, 0, 0, 0, 0, 1},
        {1, 0, 0, 0, 0, 0, 0, 0, 0, 1}, {0, 0, 0, 0, 0, 0, 1, 0, 0, 1},
        {0, 0, 0, 0, 0, 1, 0, 0, 0, 0}, {1, 0, 1, 0, 0, 0, 0, 0, 0, 0},
        {1, 0, 0, 0, 0, 0, 0, 0, 0, 1}, {0, 0, 1, 0, 0, 0, 0, 0, 0, 0},
        {0, 1, 0, 0, 0, 0, 0, 0, 0, 0}, {1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 1},
        {0, 0, 0, 1, 0, 0, 0, 0, 0, 0}, {0, 1, 0, 0, 1, 0, 1, 0, 0, 0},
        {0, 0, 1, 0, 1, 0, 0, 0, 0, 0}, {1, 0, 0, 1, 1, 0, 0, 0, 0, 1},
        {0, 1, 0, 1, 0, 0, 0, 0, 0, 0}, {0, 1, 0, 0, 0, 0, 0, 0, 0, 1},
        {0, 0, 1, 0, 1, 0, 2, 0, 0, 0}, {0, 0, 0, 0, 0, 1, 2, 0, 0, 0},
        {0, 0, 0, 0, 1, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 1, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 1, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 1, 0, 0, 0, 0}, {0, 0, 0, 1, 1, 0, 0, 0, 0, 1},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 2, 0, 0, 0, 0, 0}, {0, 2, 1, 0, 0, 0, 0, 0, 0, 1},
        {0, 0, 0, 0, 0, 1, 0, 0, 0, 0}, {0, 0, 0, 0, 1, 0, 0, 0, 0, 0},
        {1, 0, 1, 0, 1, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 1, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 1, 0, 0, 0, 0, 0, 0, 1}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 1, 1, 1, 0, 0, 1, 0, 1}, {1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 1, 0, 0, 1, 0, 0, 0, 0, 0}, {1, 0, 0, 0, 2, 0, 0, 0, 0, 2},
        {2, 0, 0, 0, 2, 0, 0, 0, 0, 0}, {1, 1, 1, 0, 0, 0, 0, 0, 0, 0},
        {0, 1, 0, 0, 0, 0, 1, 0, 0, 0}, {0, 1, 1, 0, 1, 1, 0, 0, 0, 0},
        {1, 0, 1, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 1, 1, 0, 0, 0, 0, 0},
        {0, 1, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 1, 1, 0, 0, 0, 0, 1},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {1, 0, 0, 0, 0, 1, 0, 0, 0, 1},
        {0, 0, 0, 1, 0, 0, 0, 0, 0, 0}, {1, 0, 0, 1, 0, 0, 0, 0, 0, 1},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 1, 0, 0, 0, 0, 0, 0, 0},
        {1, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 1, 0, 1, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 1}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 1, 0, 0, 0, 0, 0, 1}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 1, 0, 0, 0, 1, 0, 1, 0, 1}, {0, 0, 0, 1, 0, 0, 0, 0, 0, 1},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 1, 1, 0, 0},
        {1, 0, 0, 1, 0, 0, 1, 0, 0, 0}, {0, 1, 1, 0, 0, 0, 1, 0, 0, 1},
        {0, 0, 0, 0, 0, 1, 0, 0, 0, 1}, {0, 0, 1, 1, 1, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 1, 0, 0, 0, 0, 0}, {0, 0, 0, 1, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 1, 0, 0, 0, 0, 1}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 1, 0, 0, 0, 0, 0, 0, 1}, {0, 0, 1, 0, 1, 0, 0, 0, 0, 1},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 1}, {0, 0, 1, 0, 0, 0, 0, 0, 0, 1},
        {0, 0, 0, 0, 1, 0, 1, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 1},
        {0, 0, 0, 0, 0, 0, 0, 1, 0, 0}, {2, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {1, 0, 1, 1, 0, 1, 0, 0, 0, 1}, {0, 1, 0, 0, 0, 0, 0, 0, 0, 1},
        {0, 0, 0, 0, 0, 0, 1, 0, 0, 0}, {0, 1, 1, 0, 0, 1, 0, 0, 0, 0},
        {2, 0, 1, 0, 0, 0, 1, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 1, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 1, 1, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 1, 0, 2, 0, 0, 0, 0},
        {1, 1, 0, 0, 0, 1, 0, 0, 0, 0}, {0, 0, 1, 0, 1, 0, 0, 0, 0, 0},
        {0, 1, 0, 1, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 1, 0, 0, 0, 0, 0, 0},
        {0, 0, 1, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 1, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 1, 0, 1, 1, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 2, 0, 0, 1}, {0, 1, 0, 0, 0, 0, 0, 1, 0, 0},
        {0, 0, 0, 0, 0, 0, 1, 0, 0, 0}, {0, 0, 0, 0, 0, 1, 0, 0, 0, 1},
        {0, 1, 0, 0, 1, 0, 1, 0, 0, 0}, {0, 0, 0, 0, 0, 1, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 1, 0}, {1, 1, 0, 0, 0, 0, 0, 0, 0, 2},
        {0, 0, 1, 0, 0, 0, 0, 0, 0, 1}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 1, 0, 1, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 1, 0},
        {0, 1, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 1, 0, 0, 1, 0, 0, 1},
        {0, 0, 1, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 1, 1, 0, 0, 0}, {0, 0, 0, 0, 0, 1, 0, 0, 0, 0},
        {0, 0, 1, 0, 0, 0, 0, 0, 0, 0}, {1, 1, 0, 0, 0, 0, 0, 0, 0, 1},
        {0, 0, 0, 1, 0, 0, 0, 0, 0, 0}, {0, 1, 0, 1, 1, 0, 0, 1, 0, 0},
        {0, 0, 1, 0, 1, 0, 0, 0, 0, 1}, {0, 0, 1, 0, 0, 0, 0, 0, 0, 1},
        {0, 0, 0, 0, 1, 0, 1, 0, 0, 0}, {1, 1, 0, 0, 1, 1, 0, 0, 0, 0},
        {2, 0, 0, 1, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 1, 0, 1},
        {1, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {1, 1, 2, 0, 0, 0, 1, 0, 0, 0},
        {0, 2, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 1, 0, 0, 0, 0, 0, 0}, {0, 0, 1, 0, 0, 0, 0, 0, 0, 1},
        {0, 1, 0, 1, 1, 0, 0, 0, 0, 0}, {1, 1, 0, 0, 0, 1, 0, 0, 0, 1},
        {0, 0, 0, 1, 0, 1, 0, 0, 0, 0}, {0, 0, 0, 1, 1, 0, 0, 1, 0, 0},
        {0, 1, 0, 0, 1, 0, 2, 0, 0, 1}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 1, 0}, {1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 1}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 2},
        {1, 0, 0, 1, 1, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 1, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 1, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 1, 0, 0, 0, 1},
        {0, 0, 0, 2, 0, 0, 0, 0, 0, 0}, {0, 1, 0, 2, 1, 0, 0, 0, 0, 0},
        {0, 0, 1, 0, 0, 0, 0, 0, 0, 0}, {0, 1, 1, 0, 1, 1, 1, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 1, 0, 0, 1}, {0, 0, 0, 2, 0, 0, 0, 0, 0, 0},
        {1, 0, 1, 0, 0, 0, 0, 0, 0, 1}, {0, 0, 0, 0, 1, 0, 0, 0, 0, 1},
        {0, 0, 0, 1, 0, 0, 0, 0, 0, 0}, {0, 1, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 1, 0, 0, 2}, {0, 0, 0, 0, 2, 1, 0, 0, 0, 0},
        {0, 0, 0, 1, 1, 0, 0, 0, 0, 1}, {0, 1, 0, 1, 1, 0, 0, 0, 0, 0},
        {0, 0, 1, 0, 1, 0, 0, 0, 0, 2}, {0, 0, 1, 0, 0, 0, 1, 0, 0, 0},
        {0, 0, 0, 1, 0, 0, 0, 0, 0, 0}, {0, 0, 1, 0, 1, 0, 0, 0, 0, 0},
        {0, 1, 0, 0, 1, 0, 0, 0, 0, 1}, {1, 0, 0, 1, 1, 0, 1, 0, 0, 0},
        {0, 0, 1, 0, 0, 0, 0, 1, 0, 1}, {1, 0, 0, 0, 0, 0, 1, 0, 1, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 1, 0, 0, 0, 0, 1},
        {0, 1, 1, 1, 0, 0, 0, 0, 0, 1}, {1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {1, 0, 1, 0, 1, 1, 1, 0, 0, 1}, {0, 0, 0, 0, 2, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 1, 1, 1, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 2, 0, 0, 1, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 1, 1, 0, 2, 1, 0, 0, 0, 0}, {0, 1, 0, 0, 0, 0, 0, 0, 0, 1},
        {1, 0, 0, 1, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 1},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 1, 1, 0, 0, 0},
        {0, 0, 0, 0, 0, 1, 0, 0, 0, 0}, {0, 1, 0, 0, 1, 0, 0, 1, 0, 0},
    });
    arma::vec test_phen{
        1, 1, 1, 0, 1, 0, 0, 0, 1, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0,
        1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 0,
        0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1,
        0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 1,
        0, 1, 1, 1, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 1, 0, 0, 1, 1, 0,
        0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 0, 1, 1, 0, 0, 1,
        0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 1,
        0, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0,
        1, 1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 1, 1, 0, 0, 1, 0, 1, 0, 1,
        0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0,
        0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0,
        0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0,
        0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 1, 1, 1, 1, 0, 1, 0,
        1, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 1, 1, 0, 1,
        1, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0,
        0, 1, 0, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1,
        0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 0,
        0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0,
        0, 1, 0, 0, 0, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 0,
        0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 1, 0, 1, 0, 0,
        1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0,
        0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1,
        0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0,
        0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1,
        0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 1, 0, 1, 0, 0, 1, 1,
        0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 1, 1, 0,
        0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0,
        1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 0,
        0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1,
        1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0,
        0, 1, 1, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 0, 0,
        0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0,
        0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0,
        1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 1,
        0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0,
        1, 1, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1,
        1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0,
        1, 0, 1, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0,
        0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1,
        0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0,
        0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 1, 1, 0,
        1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 1, 1, 0};
    gene.set_matrix("test_transcript1", sim_data);
    tp.method = "CMC";
    tp.cmcmaf = 0.005;
    tp.nperm = 1000;
    tp.hotellings = true;
    cov_ptr->set_phenotype_vector(test_phen);
    {
      Methods methods(tp, cov_ptr);

      double cmc = methods.CMC(gene, cov_ptr->get_phenotype_vector(),
                               "test_transcript1", 0.005);
      REQUIRE(cmc == Approx(1.3005).epsilon(0.001));
    }

    {
      tp.nperm = 0;
      Methods methods(tp, cov_ptr);
      double cmc = methods.CMC(gene, cov_ptr->get_phenotype_vector(),
                               "test_transcript1");
      REQUIRE(cmc == Approx(0.2252).epsilon(0.001));
    }
  }
}

TEST_CASE("Gradient descent supports non-canonical links", "[glm]") {
  arma::arma_rng::set_seed(0);

  // Simple dataset with intercept and one predictor
  arma::mat X{{1.0, -2.0}, {1.0, -1.0}, {1.0, 0.0},
              {1.0, 1.0},  {1.0, 2.0},  {1.0, 3.0}};
  arma::vec Y{0, 0, 0, 1, 1, 1};

  // Fit using gradient descent with a non-canonical link
  TaskParams tp_gd;
  tp_gd.optimizer = "gradient_descent";
  Binomial cloglog_link("cloglog");
  GLM<Binomial> gd_model(X, Y, cloglog_link, tp_gd);

  REQUIRE(gd_model.success_);
  REQUIRE(arma::is_finite(gd_model.beta_));
}
