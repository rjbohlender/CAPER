#include <catch2/catch.hpp>

#include <armadillo>
#include <atomic>
#include <filesystem>
#include <limits>
#include <memory>
#include <random>
#include <sstream>
#include <string>
#include <system_error>
#include <unordered_map>
#include <vector>

#define private public
#define protected public
#include "../power/powerop.hpp"
#undef private
#undef protected

namespace {
namespace fs = std::filesystem;

std::atomic_size_t g_context_counter{0};

struct PowerOpTestContext {
  fs::path repo_root;
  fs::path output_dir;
  TaskParams tp;
  std::shared_ptr<Covariates> cov;
  Weight weight;
  std::unique_ptr<Filter> filter;
  std::unique_ptr<Gene> gene;
  std::vector<std::vector<int8_t>> permutations;
  std::unique_ptr<PowerTask> task;
  std::shared_ptr<PowerReporter> reporter;
  std::unique_ptr<PowerOp> op;

  PowerOpTestContext() {
    repo_root = fs::absolute(fs::path(__FILE__)).parent_path().parent_path();
    const auto unique_id = g_context_counter.fetch_add(1, std::memory_order_relaxed);
    output_dir = fs::temp_directory_path() /
                 fs::path("caper_powerop_sample_test-" + std::to_string(unique_id));
    fs::create_directories(output_dir);

    tp = TaskParams{};
    tp.method = "SKAT";
    tp.maf = 0.5;
    tp.mac = std::numeric_limits<arma::uword>::max();
    tp.optimizer = "irls";
    tp.no_weights = false;
    tp.aaf_filter = false;
    tp.program_directory = repo_root.string() + "/";
    tp.whitelist_path = (repo_root / "filter" / "filter_whitelist.csv").string();
    tp.output_path = output_dir.string();
    tp.power = true;
    tp.alpha = arma::vec({0.05});
    tp.bootstrap_reps = 0;
    tp.ncases = {2};
    tp.ncontrols = {2};
    tp.bin_epsilon = 0.0;

    std::stringstream test_data;
    std::stringstream test_cov;
    std::stringstream test_casm;
    std::stringstream test_ped;

    const std::string header =
        "Chr\tStart\tEnd\tRef\tAlt\tType\tGenes\tTranscripts\tRegion\tFunction\t"
        "Annotation(c.change:p.change)\tcase1\tcase2\tcontrol1\tcontrol2\n";
    test_data << header;
    test_data <<
        "chr1\t1\t1\tA\tC\tSNV\ttest_gene\ttest_transcript1\tcoding\t"
        "nonsynonymous SNV\t.\t0101\n";
    test_data <<
        "chr1\t25\t25\tA\tG\tSNV\ttest_gene\ttest_transcript1\tcoding\t"
        "nonsynonymous SNV\t.\t1111\n";
    test_data <<
        "chr1\t27\t38\tA\tATTACAGATT\tinsertion\ttest_gene\ttest_transcript1\t"
        "coding\tframeshift\t.\t0222\n";
    test_data <<
        "chr1\t55\t55\tT\tG\tSNV\ttest_gene\ttest_transcript1\texon\t"
        "nonsynonymous SNV\t.\t2001\n";
    test_data <<
        "chr1\t1\t1\tA\tC\tSNV\ttest_gene\ttest_transcript2\tcoding\t"
        "nonsynonymous SNV\t.\t0101\n";
    test_data <<
        "chr1\t25\t25\tA\tG\tSNV\ttest_gene\ttest_transcript2\tcoding\t"
        "nonsynonymous SNV\t.\t1101\n";
    test_data <<
        "chr1\t27\t38\tA\tATTACAGATT\tinsertion\ttest_gene\ttest_transcript2\t"
        "coding\tframeshift\t.\t0222\n";
    test_data <<
        "chr1\t55\t55\tT\tG\tSNV\ttest_gene\ttest_transcript2\texon\t"
        "nonsynonymous SNV\t.\t2001\n";

    test_cov << "control2\t0\t1.5\t1.5\n";
    test_cov << "case2\t1\t0.5\t0.5\n";
    test_cov << "control1\t0\t0.5\t1\n";
    test_cov << "case1\t1\t1\t0.5\n";

    test_ped << "control1\tcontrol1\t0\t0\t0\t1\n";
    test_ped << "control2\tcontrol2\t0\t0\t0\t1\n";
    test_ped << "case1\tcase1\t0\t0\t0\t2\n";
    test_ped << "case2\tcase2\t0\t0\t0\t2\n";

    test_casm <<
        "chr1\t1\t1\tA\tC\tSNV\ttest_gene\ttest_transcript1\t0.5\n";
    test_casm <<
        "chr1\t25\t25\tA\tG\tSNV\ttest_gene\ttest_transcript1\t1.5\n";
    test_casm <<
        "chr1\t27\t38\tA\tATTACAGATT\tinsertion\ttest_gene\ttest_transcript1\t3.0\n";
    test_casm <<
        "chr1\t55\t55\tT\tG\tSNV\ttest_gene\ttest_transcript1\t0.5\n";
    test_casm <<
        "chr1\t1\t1\tA\tC\tSNV\ttest_gene\ttest_transcript2\t0.5\n";
    test_casm <<
        "chr1\t25\t25\tA\tG\tSNV\ttest_gene\ttest_transcript2\t1.5\n";
    test_casm <<
        "chr1\t27\t38\tA\tATTACAGATT\tinsertion\ttest_gene\ttest_transcript2\t3.0\n";
    test_casm <<
        "chr1\t55\t55\tT\tG\tSNV\ttest_gene\ttest_transcript2\t0.5\n";

    cov = std::make_shared<Covariates>(test_ped, test_cov, tp);
    cov->sort_covariates(header);

    weight = Weight(test_casm);
    filter = std::make_unique<Filter>(
        (repo_root / "filter" / "filter_whitelist.csv").string());

    std::unordered_map<std::string, arma::uword> nvariants{
        {"test_transcript1", 4},
        {"test_transcript2", 4}};

    gene = std::make_unique<Gene>(test_data, cov, cov->get_nsamples(), nvariants,
                                  weight, tp, *filter);

    task = std::make_unique<PowerTask>(Stage::Power, *gene, cov, tp,
                                       permutations);
    reporter = std::make_shared<PowerReporter>(tp);
    op = std::make_unique<PowerOp>(*task, reporter, 1.0, false);
  }

  ~PowerOpTestContext() {
    op.reset();
    reporter.reset();
    task.reset();
    gene.reset();
    filter.reset();
    std::error_code ec;
    fs::remove_all(output_dir, ec);
  }
};
} // namespace

TEST_CASE("PowerOp::sample preserves duplicates when assembling sparse rows",
          "[powerop][sample]") {
  PowerOpTestContext ctx;
  auto &op = *ctx.op;

  arma::mat dense = {{0.0, 5.0, 0.0},
                     {1.0, 0.0, 2.0},
                     {0.0, 3.0, 4.0}};
  arma::sp_mat X(dense);

  op.cases_ = arma::uvec{1};
  op.controls_ = arma::uvec{2};
  op.gen_ = std::mt19937(17);

  arma::sp_mat sampled = op.sample(X, 2, 2);

  arma::mat expected_dense = {{1.0, 0.0, 2.0},
                              {1.0, 0.0, 2.0},
                              {0.0, 3.0, 4.0},
                              {0.0, 3.0, 4.0}};
  arma::sp_mat expected(expected_dense);

  REQUIRE(sampled.n_rows == 4);
  REQUIRE(sampled.n_cols == X.n_cols);
  REQUIRE(arma::approx_equal(arma::mat(sampled), arma::mat(expected),
                             "absdiff", 0.0));
}

TEST_CASE("PowerOp::sample copies sparse structure for random selections",
          "[powerop][sample]") {
  PowerOpTestContext ctx;
  auto &op = *ctx.op;

  arma::mat dense = {{10.0, 0.0, 1.0, 0.0},
                     {0.0, 20.0, 0.0, 1.0},
                     {2.0, 0.0, 30.0, 0.0},
                     {0.0, 0.0, 0.0, 4.0},
                     {5.0, 6.0, 0.0, 0.0}};
  arma::sp_mat X(dense);

  op.cases_ = arma::uvec{0, 2, 4};
  op.controls_ = arma::uvec{1, 3};

  const arma::uword ncases = 3;
  const arma::uword ncontrols = 2;
  constexpr unsigned int seed = 12345;
  op.gen_ = std::mt19937(seed);

  std::mt19937 reference(seed);
  std::uniform_int_distribution<> case_dis(0, op.cases_.n_elem - 1);
  std::uniform_int_distribution<> control_dis(0, op.controls_.n_elem - 1);

  std::vector<arma::uword> expected_indices;
  expected_indices.reserve(ncases + ncontrols);
  for (arma::uword i = 0; i < ncases; ++i) {
    expected_indices.push_back(op.cases_[case_dis(reference)]);
  }
  for (arma::uword i = 0; i < ncontrols; ++i) {
    expected_indices.push_back(op.controls_[control_dis(reference)]);
  }

  arma::sp_mat sampled = op.sample(X, ncases, ncontrols);

  arma::uvec expected_vec(expected_indices.size());
  for (arma::uword i = 0; i < expected_vec.n_elem; ++i) {
    expected_vec[i] = expected_indices[i];
  }

  arma::uvec actual_indices = arma::join_cols(op.cases_(op.case_idx_),
                                              op.controls_(op.control_idx_));
  REQUIRE(arma::all(actual_indices == expected_vec));

  arma::mat expected_dense(expected_vec.n_elem, X.n_cols, arma::fill::zeros);
  for (arma::uword row = 0; row < expected_vec.n_elem; ++row) {
    expected_dense.row(row) =
        arma::rowvec(arma::mat(X.row(expected_vec[row])));
  }
  arma::sp_mat expected(expected_dense);

  REQUIRE(sampled.n_rows == expected_vec.n_elem);
  REQUIRE(sampled.n_cols == X.n_cols);
  REQUIRE(arma::approx_equal(arma::mat(sampled), arma::mat(expected),
                             "absdiff", 0.0));
}
