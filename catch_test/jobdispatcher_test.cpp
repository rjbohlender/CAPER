#include <catch2/catch.hpp>
#include <fstream>
#include <filesystem>
#include <limits>
#include <memory>
#include <sstream>
#include <vector>
#include <algorithm>

#define private public
#include "../utility/jobdispatcher.hpp"
#undef private
#include "../data/filter.hpp"
#include "../caper/capertask.hpp" // for Stage

class DummyTask {
public:
    std::string gene;
    arma::uword success;
    arma::uword perm_count;
    int perm_offset;
    DummyTask(Stage, Gene &g, std::shared_ptr<Covariates>, TaskParams &, std::vector<std::vector<int8_t>> &)
        : gene(g.gene_name), success(0), perm_count(0), perm_offset(0) {}
    DummyTask(Stage, Gene &g, std::shared_ptr<Covariates>, TaskParams &, arma::uword s, arma::uword pc, int off, arma::uword, std::vector<std::vector<int8_t>> &)
        : gene(g.gene_name), success(s), perm_count(pc), perm_offset(off) {}
};

class DummyReporter {
public:
    std::vector<std::string> genes;
    std::vector<DummyTask> tasks;
    void set_ncases(arma::uword) {}
    void set_ncontrols(arma::uword) {}
    template <typename T>
    void report(const T &, const TaskParams &) {}
    template <typename T>
    void vaast_sample_index_map(const T &) {}
    void cleanup(const TaskParams &) {}
};

class DummyOp {
public:
    bool done_;
    DummyTask carvaTask;
    std::shared_ptr<DummyReporter> reporter_;

    DummyOp(DummyTask &t, std::shared_ptr<DummyReporter> reporter, double, bool)
        : done_(true), carvaTask(t), reporter_(std::move(reporter)) {
        reporter_->genes.push_back(carvaTask.gene);
        reporter_->tasks.push_back(carvaTask);
    }
    DummyOp(DummyTask &&t, std::shared_ptr<DummyReporter> reporter, double, bool)
        : done_(true), carvaTask(std::move(t)), reporter_(std::move(reporter)) {
        reporter_->genes.push_back(carvaTask.gene);
        reporter_->tasks.push_back(carvaTask);
    }
    void run() {}
    void finish() {}
};

TEST_CASE("JobDispatcher dispatches tasks for each gene") {
    namespace fs = std::filesystem;
    auto tmp = fs::temp_directory_path();

    auto ped_path = tmp / "jd_ped.ped";
    std::ofstream ped(ped_path);
    ped << "#FID\tIID\tFather\tMother\tSex\tPhenotype\n";
    ped << "control1\tcontrol1\t0\t0\t0\t1\n";
    ped << "control2\tcontrol2\t0\t0\t0\t1\n";
    ped << "case1\tcase1\t0\t0\t0\t2\n";
    ped << "case2\tcase2\t0\t0\t0\t2\n";
    ped.close();

    auto input_path = tmp / "jd_input.tsv";
    std::ofstream input(input_path);
    input << "Chr\tStart\tEnd\tRef\tAlt\tType\tGenes\tTranscripts\tRegion\tFunction\tAnnotation(c.change:p.change)\tcase1\tcase2\tcontrol1\tcontrol2\n";
    input.close();

    TaskParams tp{};
    tp.covariates_path.clear();
    tp.ped_path = ped_path.string();
    tp.input_path = input_path.string();
    tp.whitelist_path = (fs::path(__FILE__).parent_path().parent_path() / "filter" / "filter_whitelist.csv").string();
    tp.nthreads = 2;
    tp.nperm = 0;
    tp.max_perms = 0; // prevent constructor from auto-dispatching
    tp.mac = std::numeric_limits<arma::uword>::max();
    tp.maf = 1.0;
    tp.min_variant_count = 0;
    tp.min_minor_allele_count = 0;
    tp.no_weights = true;
    tp.nocovadj = true;
    tp.optimizer = "irls";
    tp.method = "BURDEN";

    auto reporter = std::make_shared<DummyReporter>();
    JobDispatcher<DummyOp, DummyTask, DummyReporter> jd(tp, reporter);

    // JobDispatcher frees its covariates after construction; reinitialize for manual dispatch
    jd.cov_ = std::make_shared<Covariates>(tp);
    jd.cov_->sort_covariates(jd.header_);

    std::stringstream ss;
    ss << "chr1\t1\t1\tA\tG\tSNV\tGene1\tTranscript1\tcoding\tnonsynonymous SNV\t.\t0101\n";
    ss << "chr1\t2\t2\tT\tC\tSNV\tGene2\tTranscript1\tcoding\tnonsynonymous SNV\t.\t0101\n";

    Filter filter(tp.whitelist_path);
    jd.all_gene_dispatcher(ss, filter);

    REQUIRE(jd.tq_.size() == 2);
    REQUIRE(reporter->genes.size() == 2);
    REQUIRE(reporter->genes[0] == "Gene1");
    REQUIRE(reporter->genes[1] == "Gene2");
}

TEST_CASE("JobDispatcher loads external permutations from file") {
    namespace fs = std::filesystem;
    auto tmp = fs::temp_directory_path();

    auto ped_path = tmp / "jd_ped.ped";
    std::ofstream ped(ped_path);
    ped << "#FID\tIID\tFather\tMother\tSex\tPhenotype\n";
    ped << "control1\tcontrol1\t0\t0\t0\t1\n";
    ped << "control2\tcontrol2\t0\t0\t0\t1\n";
    ped << "case1\tcase1\t0\t0\t0\t2\n";
    ped << "case2\tcase2\t0\t0\t0\t2\n";
    ped.close();

    auto input_path = tmp / "jd_input.tsv";
    std::ofstream input(input_path);
    input << "Chr\tStart\tEnd\tRef\tAlt\tType\tGenes\tTranscripts\tRegion\tFunction\tAnnotation(c.change:p.change)\tcase1\tcase2\tcontrol1\tcontrol2\n";
    input.close();

    auto external_path = tmp / "jd_external_perms.txt";
    std::ofstream external(external_path);
    external << "0101\n";
    external << "1010\n";
    external.close();

    TaskParams tp{};
    tp.covariates_path.clear();
    tp.ped_path = ped_path.string();
    tp.input_path = input_path.string();
    tp.whitelist_path = (fs::path(__FILE__).parent_path().parent_path() / "filter" / "filter_whitelist.csv").string();
    tp.nthreads = 2;
    tp.nperm = 0;
    tp.max_perms.reset();
    tp.mac = std::numeric_limits<arma::uword>::max();
    tp.maf = 1.0;
    tp.min_variant_count = 0;
    tp.min_minor_allele_count = 0;
    tp.no_weights = true;
    tp.nocovadj = true;
    tp.optimizer = "irls";
    tp.method = "BURDEN";
    tp.external = true;
    tp.external_path = external_path.string();

    auto reporter = std::make_shared<DummyReporter>();
    JobDispatcher<DummyOp, DummyTask, DummyReporter> jd(tp, reporter);

    REQUIRE(jd.permutation_ptr_ != nullptr);
    REQUIRE(jd.permutation_ptr_->size() == 2);
    REQUIRE(jd.permutation_ptr_->at(0) == std::vector<int8_t>{0, 1, 0, 1});
    REQUIRE(jd.permutation_ptr_->at(1) == std::vector<int8_t>{1, 0, 1, 0});
    REQUIRE(jd.tp_.nperm == 2);
}

TEST_CASE("gene_list_dispatcher only dispatches listed genes") {
    namespace fs = std::filesystem;
    auto tmp = fs::temp_directory_path();

    auto ped_path = tmp / "jd_ped.ped";
    std::ofstream ped(ped_path);
    ped << "#FID\tIID\tFather\tMother\tSex\tPhenotype\n";
    ped << "control1\tcontrol1\t0\t0\t0\t1\n";
    ped << "control2\tcontrol2\t0\t0\t0\t1\n";
    ped << "case1\tcase1\t0\t0\t0\t2\n";
    ped << "case2\tcase2\t0\t0\t0\t2\n";
    ped.close();

    auto input_path = tmp / "jd_input.tsv";
    std::ofstream input(input_path);
    input << "Chr\tStart\tEnd\tRef\tAlt\tType\tGenes\tTranscripts\tRegion\tFunction\tAnnotation(c.change:p.change)\tcase1\tcase2\tcontrol1\tcontrol2\n";
    input.close();

    TaskParams tp{};
    tp.covariates_path.clear();
    tp.ped_path = ped_path.string();
    tp.input_path = input_path.string();
    tp.whitelist_path = (fs::path(__FILE__).parent_path().parent_path() / "filter" / "filter_whitelist.csv").string();
    tp.nthreads = 2;
    tp.nperm = 0;
    tp.max_perms = 0; // prevent constructor from auto-dispatching
    tp.mac = std::numeric_limits<arma::uword>::max();
    tp.maf = 1.0;
    tp.min_variant_count = 0;
    tp.min_minor_allele_count = 0;
    tp.no_weights = true;
    tp.nocovadj = true;
    tp.optimizer = "irls";
    tp.method = "BURDEN";
    tp.gene_list = std::string("Gene2");

    auto reporter = std::make_shared<DummyReporter>();
    JobDispatcher<DummyOp, DummyTask, DummyReporter> jd(tp, reporter);

    // JobDispatcher frees its covariates after construction; reinitialize for manual dispatch
    jd.cov_ = std::make_shared<Covariates>(tp);
    jd.cov_->sort_covariates(jd.header_);

    std::stringstream ss;
    ss << "chr1\t1\t1\tA\tG\tSNV\tGene1\tTranscript1\tcoding\tnonsynonymous SNV\t.\t0101\n";
    ss << "chr1\t2\t2\tT\tC\tSNV\tGene2\tTranscript2\tcoding\tnonsynonymous SNV\t.\t0101\n";
    ss << "chr1\t3\t3\tG\tA\tSNV\tGene3\tTranscript3\tcoding\tnonsynonymous SNV\t.\t0101\n";

    Filter filter(tp.whitelist_path);
    jd.gene_list_dispatcher(ss, filter);

    REQUIRE(jd.tq_.size() == 1);
    REQUIRE(jd.ngenes_ == 1);
    REQUIRE(jd.gene_list_.empty());
    REQUIRE(reporter->tasks.size() == 1);
    REQUIRE(reporter->genes.size() == 1);
    REQUIRE(reporter->genes[0] == "Gene2");
}

TEST_CASE("JobDispatcher throws when genes are not sorted") {
    namespace fs = std::filesystem;
    auto tmp = fs::temp_directory_path();

    auto ped_path = tmp / "jd_ped.ped";
    std::ofstream ped(ped_path);
    ped << "#FID\tIID\tFather\tMother\tSex\tPhenotype\n";
    ped << "control1\tcontrol1\t0\t0\t0\t1\n";
    ped << "control2\tcontrol2\t0\t0\t0\t1\n";
    ped << "case1\tcase1\t0\t0\t0\t2\n";
    ped << "case2\tcase2\t0\t0\t0\t2\n";
    ped.close();

    auto input_path = tmp / "jd_input.tsv";
    std::ofstream input(input_path);
    input << "Chr\tStart\tEnd\tRef\tAlt\tType\tGenes\tTranscripts\tRegion\tFunction\tAnnotation(c.change:p.change)\tcase1\tcase2\tcontrol1\tcontrol2\n";
    input.close();

    TaskParams tp{};
    tp.covariates_path.clear();
    tp.ped_path = ped_path.string();
    tp.input_path = input_path.string();
    tp.whitelist_path = (fs::path(__FILE__).parent_path().parent_path() / "filter" / "filter_whitelist.csv").string();
    tp.nthreads = 2;
    tp.nperm = 0;
    tp.max_perms = 0; // prevent constructor from auto-dispatching
    tp.mac = std::numeric_limits<arma::uword>::max();
    tp.maf = 1.0;
    tp.min_variant_count = 0;
    tp.min_minor_allele_count = 0;
    tp.no_weights = true;
    tp.nocovadj = true;
    tp.optimizer = "irls";
    tp.method = "BURDEN";

    auto reporter = std::make_shared<DummyReporter>();
    JobDispatcher<DummyOp, DummyTask, DummyReporter> jd(tp, reporter);

    // JobDispatcher frees its covariates after construction; reinitialize for manual dispatch
    jd.cov_ = std::make_shared<Covariates>(tp);
    jd.cov_->sort_covariates(jd.header_);

    std::stringstream ss;
    ss << "chr1\t1\t1\tA\tG\tSNV\tGene1\tTranscript1\tcoding\tnonsynonymous SNV\t.\t0101\n";
    ss << "chr1\t2\t2\tT\tC\tSNV\tGene2\tTranscript1\tcoding\tnonsynonymous SNV\t.\t0101\n";
    ss << "chr1\t3\t3\tG\tA\tSNV\tGene1\tTranscript1\tcoding\tnonsynonymous SNV\t.\t0101\n";

    Filter filter(tp.whitelist_path);
    REQUIRE_THROWS_WITH(
        jd.all_gene_dispatcher(ss, filter),
        "Gene list must be sorted by gene name. Gene Gene1 appears again on line 2 of the gene stream. Please sort the gene stream by gene name and transcript.");
}

TEST_CASE("JobDispatcher skips genes when variants are masked by BED") {
    namespace fs = std::filesystem;
    auto tmp = fs::temp_directory_path();

    auto ped_path = tmp / "jd_ped.ped";
    std::ofstream ped(ped_path);
    ped << "#FID\tIID\tFather\tMother\tSex\tPhenotype\n";
    ped << "control1\tcontrol1\t0\t0\t0\t1\n";
    ped << "control2\tcontrol2\t0\t0\t0\t1\n";
    ped << "case1\tcase1\t0\t0\t0\t2\n";
    ped << "case2\tcase2\t0\t0\t0\t2\n";
    ped.close();

    auto input_path = tmp / "jd_input.tsv";
    std::ofstream input(input_path);
    input << "Chr\tStart\tEnd\tRef\tAlt\tType\tGenes\tTranscripts\tRegion\tFunction\tAnnotation(c.change:p.change)\tcase1\tcase2\tcontrol1\tcontrol2\n";
    input.close();

    auto bed_path = tmp / "jd_mask.bed";
    std::ofstream bed(bed_path);
    bed << "chr1\t1\t1\tA\tG\n";
    bed << "chr1\t2\t2\tT\tC\n";
    bed.close();

    TaskParams tp{};
    tp.covariates_path.clear();
    tp.ped_path = ped_path.string();
    tp.input_path = input_path.string();
    tp.whitelist_path = (fs::path(__FILE__).parent_path().parent_path() / "filter" / "filter_whitelist.csv").string();
    tp.nthreads = 2;
    tp.nperm = 0;
    tp.max_perms = 0; // prevent constructor from auto-dispatching
    tp.mac = std::numeric_limits<arma::uword>::max();
    tp.maf = 1.0;
    tp.min_variant_count = 0;
    tp.min_minor_allele_count = 0;
    tp.no_weights = true;
    tp.nocovadj = true;
    tp.optimizer = "irls";
    tp.method = "BURDEN";
    tp.bed = bed_path.string();

    auto reporter = std::make_shared<DummyReporter>();
    JobDispatcher<DummyOp, DummyTask, DummyReporter> jd(tp, reporter);

    // JobDispatcher frees its covariates after construction; reinitialize for manual dispatch
    jd.cov_ = std::make_shared<Covariates>(tp);
    jd.cov_->sort_covariates(jd.header_);

    std::stringstream ss;
    ss << "chr1\t1\t1\tA\tG\tSNV\tGene1\tTranscript1\tcoding\tnonsynonymous SNV\t.\t0101\n";
    ss << "chr1\t2\t2\tT\tC\tSNV\tGene1\tTranscript1\tcoding\tnonsynonymous SNV\t.\t0101\n";

    Filter filter(tp.whitelist_path);
    jd.all_gene_dispatcher(ss, filter);

    REQUIRE(jd.tq_.size() == 0);
    REQUIRE(jd.ngenes_ == 0);
    REQUIRE(reporter->genes.empty());
}

TEST_CASE("JobDispatcher counts only unmasked variants with BED mask") {
    namespace fs = std::filesystem;
    auto tmp = fs::temp_directory_path();

    auto ped_path = tmp / "jd_ped.ped";
    std::ofstream ped(ped_path);
    ped << "#FID\tIID\tFather\tMother\tSex\tPhenotype\n";
    ped << "control1\tcontrol1\t0\t0\t0\t1\n";
    ped << "control2\tcontrol2\t0\t0\t0\t1\n";
    ped << "case1\tcase1\t0\t0\t0\t2\n";
    ped << "case2\tcase2\t0\t0\t0\t2\n";
    ped.close();

    auto input_path = tmp / "jd_input.tsv";
    std::ofstream input(input_path);
    input << "Chr\tStart\tEnd\tRef\tAlt\tType\tGenes\tTranscripts\tRegion\tFunction\tAnnotation(c.change:p.change)\tcase1\tcase2\tcontrol1\tcontrol2\n";
    input.close();

    auto bed_path = tmp / "jd_partial_mask.bed";
    std::ofstream bed(bed_path);
    bed << "chr1\t1\t1\tA\tG\n";
    bed.close();

    TaskParams tp{};
    tp.covariates_path.clear();
    tp.ped_path = ped_path.string();
    tp.input_path = input_path.string();
    tp.whitelist_path = (fs::path(__FILE__).parent_path().parent_path() / "filter" / "filter_whitelist.csv").string();
    tp.nthreads = 2;
    tp.nperm = 0;
    tp.max_perms = 0; // prevent constructor from auto-dispatching
    tp.mac = std::numeric_limits<arma::uword>::max();
    tp.maf = 1.0;
    tp.min_variant_count = 0;
    tp.min_minor_allele_count = 0;
    tp.no_weights = true;
    tp.nocovadj = true;
    tp.optimizer = "irls";
    tp.method = "BURDEN";
    tp.bed = bed_path.string();

    auto reporter = std::make_shared<DummyReporter>();
    JobDispatcher<DummyOp, DummyTask, DummyReporter> jd(tp, reporter);

    // JobDispatcher frees its covariates after construction; reinitialize for manual dispatch
    jd.cov_ = std::make_shared<Covariates>(tp);
    jd.cov_->sort_covariates(jd.header_);

    std::stringstream ss;
    ss << "chr1\t1\t1\tA\tG\tSNV\tGene1\tTranscript1\tcoding\tnonsynonymous SNV\t.\t0101\n";
    ss << "chr1\t2\t2\tT\tC\tSNV\tGene1\tTranscript1\tcoding\tnonsynonymous SNV\t.\t0101\n";

    Filter filter(tp.whitelist_path);
    jd.all_gene_dispatcher(ss, filter);

    REQUIRE(jd.tq_.size() == 1);
    REQUIRE(jd.nvariants_.count("Transcript1") == 1);
    REQUIRE(jd.nvariants_.at("Transcript1") == 1);
    REQUIRE(jd.nvariants_.size() == 1);
    REQUIRE(jd.ngenes_ == 1);
    REQUIRE(reporter->tasks.size() == 1);
    REQUIRE(reporter->genes.size() == 1);
    REQUIRE(reporter->genes[0] == "Gene1");
}

TEST_CASE("multiple_dispatch splits permutation ranges without overlap") {
    namespace fs = std::filesystem;
    auto tmp = fs::temp_directory_path();

    auto ped_path = tmp / "jd_ped.ped";
    std::ofstream ped(ped_path);
    ped << "#FID\tIID\tFather\tMother\tSex\tPhenotype\n";
    ped << "control1\tcontrol1\t0\t0\t0\t1\n";
    ped << "control2\tcontrol2\t0\t0\t0\t1\n";
    ped << "case1\tcase1\t0\t0\t0\t2\n";
    ped << "case2\tcase2\t0\t0\t0\t2\n";
    ped.close();

    auto input_path = tmp / "jd_input.tsv";
    std::ofstream input(input_path);
    input << "Chr\tStart\tEnd\tRef\tAlt\tType\tGenes\tTranscripts\tRegion\tFunction\tAnnotation(c.change:p.change)\tcase1\tcase2\tcontrol1\tcontrol2\n";
    input.close();

    TaskParams tp{};
    tp.covariates_path.clear();
    tp.ped_path = ped_path.string();
    tp.input_path = input_path.string();
    tp.whitelist_path = (fs::path(__FILE__).parent_path().parent_path() / "filter" / "filter_whitelist.csv").string();
    tp.nthreads = 3;
    tp.nperm = 11;
    tp.success_threshold = 7;
    tp.mac = std::numeric_limits<arma::uword>::max();
    tp.maf = 1.0;
    tp.min_variant_count = 0;
    tp.min_minor_allele_count = 0;
    tp.no_weights = true;
    tp.nocovadj = true;
    tp.optimizer = "irls";
    tp.method = "BURDEN";

    auto reporter = std::make_shared<DummyReporter>();
    JobDispatcher<DummyOp, DummyTask, DummyReporter> jd(tp, reporter);

    jd.cov_ = std::make_shared<Covariates>(tp);
    jd.cov_->sort_covariates(jd.header_);

    std::stringstream gene_stream;
    gene_stream << jd.header_ << "\n";
    gene_stream << "chr1\t1\t1\tA\tG\tSNV\tGeneX\tTranscript1\tcoding\tnonsynonymous SNV\t.\t0101\n";

    Filter filter(tp.whitelist_path);
    jd.nvariants_.clear();
    jd.nvariants_["Transcript1"] = 1;
    Gene gene(gene_stream, jd.cov_, jd.cov_->get_nsamples(), jd.nvariants_, jd.weight_, tp, filter);

    jd.multiple_dispatch(gene);

    REQUIRE(reporter->tasks.size() == jd.tq_.get_nthreads());

    std::vector<std::pair<int, arma::uword>> ranges;
    arma::uword total_perm = 0;
    arma::uword total_success = 0;
    for (const auto &t : reporter->tasks) {
        ranges.emplace_back(t.perm_offset, t.perm_offset + t.perm_count);
        total_perm += t.perm_count;
        total_success += t.success;
    }
    std::sort(ranges.begin(), ranges.end());
    for (size_t i = 1; i < ranges.size(); ++i) {
        REQUIRE(ranges[i - 1].second <= ranges[i].first);
    }
    REQUIRE(total_perm == tp.nperm);
    REQUIRE(total_success == tp.success_threshold);
}

