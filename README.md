# Covariate Adjusted PERmutation (CAPER) #

CAPER is a tool for gene-based testing of rare variants in sequencing studies of any size. We provide methods for the
analysis of data provided in simple, plain text, array formats. Permutation and statistical computation is
multithreaded.

This software has been carefully designed to reduce memory overhead as much as possible through the use of shared read
only objects, and elimination of data that is no longer needed. Total concurrent memory usage increases with the number
of threads. By default, the software uses half of the available number of cpus, but this can be controlled.

All methods are available for the analysis of binary traits. Quantitative traits can be analyzed using the BURDEN, RVT1,
RVT2, SKAT, and SKATO methods currently.

## Approach ##

We have developed a simple and fast approach to approximate covariate adjusted permutation. In short, we use covariate
adjusted permutation while binning samples by the similarity of their odds of being a case. We permute the phenotype of
all individuals according to their odds estimated from logistic regression. All methods can be used with covariate
adjusted permutation. If the user doesn't supply a covariate file, or disables covariate adjustment, the phenotypes are
shuffled uniformly.

## Supported Methods ##

Methods can be chosen using the -m, or --method option. The default method is VAAST.

	- BURDEN (Wu, Guan, Pankow 2016)
	- CALPHA (Neale et al. 2011)
	- CMC (Li, Leal 2008)
	- CMC1df -- OR as test statistic
	- RVT1 (Morris, Zeggini 2010)
	- RVT2 (Morris, Zeggini 2010)
	- SKAT (Wu et al. 2011; Wu, Guan, Pankow 2016)
	- SKATO (Lee, Wu, Lin 2012; Wu, Guan, Pankow 2016)
	- VAAST (Yandell et al. 2011) -- Default
	- VT (Price et al. 2010)
	- WSS (Madsen, Browning 2009)

A subset of methods can provide analytic p-values if run with --nperm 0. Those include, BURDEN, CMC, RVT1, RVT2, SKAT,
and SKATO.

### SKAT / SKAT-O ###

SKAT and SKATO are implemented using the method described by Wu, Guan, and Pankow (2016). The variant based statistic
provides significant computational speedup over the individual based statistic for very large datasets. Additionally,
this allows SKAT to be used on very large datasets, where otherwise the NxN covariance matrix will exhaust available
memory. Even with the computation advantage that this approach provides, it is still time-consuming to permute SKAT-O.

Note that SKAT can return two different values. When used without permutation, SKAT returns an analytic p-value. When
used with permutation, the costly calculation of the p-value is skipped, and the test statistic is returned instead.

## Compiling ##

Dependencies: 
- Tested with Armadillo >= 8.600
- C++ compiler supporting C++20
- C++ Boost Library > 1.66 (required for quadrature)

The dependencies can be installed via your package manager. Using Homebrew on macOS:

```bash
brew install armadillo
brew install boost

```

If you're working in a cluster or otherwise managed environment, ensure that Armadillo is compiled against both Lapack
and BLAS, or an alternative like Intel MKL. Without that, this software will fail when calling Singular Value
Decomposition and other algorithms not provided by Armadillo itself.

Create a build directory and run cmake.

```bash
mkdir build
cd build

cmake -DCMAKE_BUILD_TYPE=Release ..
make
```

If cmake fails to detect Armadillo, but you're sure it is available, you may need to direct cmake to the library, e.g.,
when compiling on a cluster, with packages in non-standard locations. In that case the following should work:

```bash
mkdir build && cd build

cmake -DCMAKE_BUILD_TYPE=Release -DARMADILLO_INCLUDE_DIR=<path_to_armadillo>/include/ -DARMADILLO_LIBRARY=<path_to_armadillo>/lib64/libarmadillo.so
make
```

If you're having trouble with cmake detecting the correct compiler.

```bash
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=<path_to_compiler> ..
make

```

The location for boost may need to be specified if it isn't installed in a typical location.

```bash
cmake -DBOOST_ROOT=<path_to_boost> ..
```

You can combine the above as necessary. Earlier versions of the Armadillo library may work, but haven't been tested. If
you need to change the compiler used from the one automatically detected to another, perhaps newer compiler:

```bash
cmake -DCMAKE_CXX_COMPILER=<path_to_executable> ..
```

## Benchmarks ##

Testing on a mid-2015 MacBook Pro, 2.2Ghz Core i7, 16GB 1600 Mhz DDR3 ram.

|   Test | Samples | Genes | Permutations | Time (sec) |
|-------:|:-------:|:-----:|:------------:|:----------:|
|    CMC | 100000  | 1000  |      0       |   341.42   |
|   RVT1 | 100000  | 1000  |      0       |   377.89   |
|   RVT2 | 100000  | 1000  |      0       |   405.56   |
|  VAAST | 100000  | 1000  |    10000     |  2086.53   |
| SKAT-O | 100000  | 1000  |      0       |   487.79   |
|    WSS | 100000  | 1000  |    10000     |  3432.29   |
|   SKAT | 100000  | 1000  |    10000     |  28293.36  |

### Gzipped Data ###

1049.7s with gzipped data, 1000 permutations of VAAST.
1063.4s with unzipped data, 1000 permutations of VAAST.

Working with gzipped data has no significant impact on runtime.

## File Formats ##

The file formats used are simple, plain text formats, which may be gzipped.

### Matrix ###

Matrix files may be zipped with gzip, or unzipped. They are provided to the program with the "-i" option. The matrix
format used for the genotype file includes a header and is as follows:

	1) Chromosome
	2) Start position (bp)
	3) End position (bp)
	4) Type (SNV) 
	5) Reference allele
	6) Alternate allele
	7) Gene symbol (e.g., BRCA1)
	8) Transcript (e.g., NM_700030) 
	9) Region (e.g., exonic, intronic)
	10) Function (synonymous, nonsynonymous)
	11+) Sample genotype (e.g., 0/1/2 - assumes alternate allele count)

Note that the type field is used for grouping in VAAST, and the program will fail if the annotation is not present when
grouping is used.

Variants are expected to be repeated for each transcript they appear in. This does increase the size of the file, but
simplifies parsing. Because the format is simple, and QC is expected to be finished before this program is run, the file
remains reasonably lean. E.g., with 20,000 simulated genes, and 100,000 samples, a matrix file, uncompressed, is only
56GB. Compressed with gzip defaults, the same matrix file is only 841MB.

### Bed Mask File ###

The "-b" option allows the user to provide a ".bed" format file to mask problematic variants. The .bed format for the
mask file is as follows:

	1) Chromosome
	2) Start position of masked region
	3) End position of masked region
    4) Reference allele
    5) Alternate allele

XQC (in the Cross Platform Association Toolkit (XPAT)) will provide a mask file for those variants that fail QC if it is
being used.

### Covariates ###

Covariates can be provided in a matrix format file. The file format is as follows:

    1) Sample ID
    2) First covariate
    3) Second covariate
    4) ...

There is no limit on the number of covariates that can be provided. If a given covariate cannot be converted to a
floating point representation, then it is assumed to be a discrete category (e.g., male/female), and will be separated
into n-1 0/1 variables where n is the number of categories.

### Weights ###

The file format for weights. The columns are required to uniquely identify variants in the case of weights for variants
with different annotations (e.g., a splice variant in one transcript and missense in another). Weights will be used by
VAAST, and SKAT / SKAT-O if the linear kernel is used.

    1) Chromosome
    2) Start position
    3) End position
    4) Reference Allele
    5) Alternate Allele
    6) Type (e.g. SNV, Deletion, Insertion, etc.)
    7) Gene symbol
    8) Transcript
    9) Weight

Weights provided via this option overwrite the beta distribution values, which are used for weighting by default in
SKAT / SKAT-O. Passing the --no_weights option will remove all weights. Using SKAT or SKAT-O with the weighted Linear
kernel option will overwrite the weights with those calculated from a beta distribution.

## Filter Whitelist ##

This is a matrix of variant types and functions, to include for each method. Any types or functions not in the list will
be ignored / removed during parsing. The whitelist format is simple. A default whitelist is provided in the filter
directory.

The default allows all variant types and functions for all methods. Removing the SPDA line would result in all variants
labelled as splicing donor / acceptor being removed.

## Testability ##

We compare the score of the most extreme phenotype distribution against the distribution of permuted statistics for
methods without an analytical p-value. The testability of this gene is then determined based upon the achievable p-value
for the gene.

## Usage ##

The program is run from the command line. The following options are available:

```
Permutation tool for gene-based rare-variant analysis.
Allowed options:

Required:
    -i [ --input ] arg                    Genotype matrix file path.
    -p [ --ped ] arg                      Path to the .ped file containing the sample phenotypes.
    -o [ --output ] arg                   Path to output directory. Two files will be output: a simple transcript level results file, and a detailed variant level result file.

Optional:
    -c [ --covariates ] arg               The covariate matrix file, tab or space separated. Format = sample_id cov1 ...
    -b [ --bed_filter ] arg               A bed file, or a comma separated list of bed files, to be used as a filter. All specified variants will be excluded.
    -f [ --filter ] arg                   A csv whitelist of TYPE and FUNCTION annotations. Default whitelist can be found in the filter directory.
    -w [ --weights ] arg                  A file providing weights.
    --no_weights                          Disable weights. Disables default CASM weights.
    --impute_to_mean                      Impute the to mean AF of cases for case samples, and the mean AF of controls for control samples.
    --whole_gene                          Analyze the union of all transcripts for a gene.
    -t [ --nthreads ]                     The number of threads. Minimum number of threads = 2. n + 1 threads, with one parent thread and n threads processing genes.
    -m [ --method ] arg (=VAAST)          The statistical method to be used. Options: {BURDEN, CALPHA, CMC, CMC1df, RVT1, RVT2, SKAT, SKATO, SKATC, VAAST, VT, WSS}.
    --optimizer arg (=irls)               The optimizer used to fit the GLM. Options: {irls, irls_svdnewton, irls_qr, irls_qr_R, gradient_descent}.
    --range arg                           A range of genes to analyze from the matrix file. Takes two values, a start gene number, and end gene number. The program will only provide results for the values in that range. Gene count starts at 1. Useful for starting multiple jobs on a cluster each processing part of a file. Note: Somewhat slower than splitting the input matrix file.
    --nperm arg (=10000)                  The maximum number of permutations to be performed.
    --mac arg                             Alternative or minor allele count cutoff per variant.
    -r [ --maf ] arg (=0.5)               Alternative or minor allele frequency cutoff per variant. We recommend using an external sample and filtering variants based on the frequency in that sample, rather than filtering within. Can result in a reduction in power for variants near the threshold.
    --aaf                                 Filter variants where the alternate allele is the major allele.
    -j [ --pthresh ] arg                  The threshold to terminate permutation based on whether it is outside the p-value CI.
    --top_only                            Output only the top transcript in the simple file.
    -s [ --successes ] arg (=200)         Number of successes for early termination.
    -l [ --genes ] arg                    A comma-separated list of genes to analyze.
    --no_detail                           Don't produce detailed, variant level output.
    --output_stats                        Write permuted statistics to .simple file following default output.
    --permute_out arg                     Output permutations to the given file. Exits after generating permutations.
    --min_minor_allele_count arg (=1)     Minimum number of minor allele copies to test a gene.
    --min_variant_count arg (=1)          Minimum number of variants to test a gene.
    --max_levels arg (=100)               Maximum number of levels for a single variable. Will be split into n-1 dummy variables.
    --bin_epsilon arg (=0.0001)           Odds closer together than the given value will be collapsed into a single bin for permutation.
    --max_perms arg                       Maximum number of permutations, used in combination with --nperm to manage memory usage. Run permutation in blocks of size nperm, up to the maximum set here. Only genes requiring additional permutation will be permuted. If you are running a small number of permutations, do not set this option.
    --seed arg                            A defined seed passed to the random number generators used for each gene.
    --check_testability arg               Return results for genes with a minimum achievable p-value less than or equal to what is given.

VAAST Options:
    -g [ --group_size ] arg (=4)          Group size, minor allele count threshold for grouping a variant. VAAST can collapse variants into groups of variants, dependent upon the collapse having a higher total VAAST score.
    --soft_maf_filter arg                 Caps the highest allele frequency for the control set in the likelihood calculation. Penalizes common variants without removing them.
    --biallelic                           Additional term for biallelic variants. For detecting potentially recessive variants.
    --site_penalty arg (=2)               VAAST site penalty. AIC penalty applied to each site in VAAST.
    --alternate_grouping                  If enabled variants are grouped all together, otherwise by VAAST 2.0 type annotation.

SKAT Options:
    -k [ --kernel ] arg (=wLinear)        Kernel for use with SKAT. One of: {wLinear, Linear}.
    --qtl                                 Analyze a quantitative trait. Values are assumed to be finite floating point values.
    --beta_weights arg (=1,25)            Parameters for the beta distribution. Two values, comma separated corresponding to a,b.
    --saddlepoint                         Force the saddlepoint approximation. Useful for highly skewed case/control sample sizes.

CMC Options:
    --cmcmaf arg (=0.005)                 Minor allele frequency cutoff for CMC collapsing.
    --hotellings                          Use Hotellings T2 instead of a chi-square test.

```

Here are some examples of how to run the program:

A simple run with 10,000 permutations and no covariates.
```bash
./caper -i matrix.txt -p phenotypes.ped -o output_dir -m VAAST -n 10000
```

A run with covariates and 1,000,000 permutations. This run will generate permutations in sets of 10,000 to prevent
exhaustion of memory.
```bash
./caper -i matrix.txt -p phenotypes.ped -o output_dir -m VAAST -n 10000 --max_perms 1000000 -c covariates.txt
```

A run with a gene range, and a seed for the random number generator. This will only analyze genes 100 to 200 in the matrix file.
```bash
./caper -i matrix.txt -p phenotypes.ped -o output_dir -m VAAST -n 10000 --range 100 200 --seed 12345
```

A run with a gene list. This will only analyze the genes in the list.
```bash
./caper -i matrix.txt -p phenotypes.ped -o output_dir -m SKAT -n 10000 -l gene1,gene2,gene3
```

A run with a limited bin epsilon. This will collapse bins with odds closer than 0.00001 together, increasing the speed of permutation.
```bash
./caper -i matrix.txt -p phenotypes.ped -o output_dir -m RVT1 -n 10000 --bin_epsilon 0.00001
```