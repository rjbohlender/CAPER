# Covariate Adjusted Rare Variant Association (CARVA) #

CARVA is a tool for gene-based testing of rare variants in
sequencing studies of any size. We provide methods for the analysis of
data provided in simple, plain text, array formats. Permutation and
statistical computation is multithreaded.

This software has been carefully designed to reduce memory overhead as
much as possible through the use of shared read only objects, and
elimination of data that is no longer needed. Total concurrent memory
usage increases with the number of threads. By default, the software
uses half of the available number of cpus, but this can be controlled.

All methods are available for the analysis of binary
traits. Quantitative traits can be analyzed using the BURDEN, RVT1,
RVT2, SKAT, and SKATO methods currently.

## Approach ##

We have implemented an alternative to the efficient resampling
approach presented in Lee, Fuchsberger, Kim, and Scott (2016). In
short, we use covariate adjusted permutation while binning samples by the
similarity of their odds of being a case.  We permute the phenotype of all 
individuals according to their odds estimated from logistic regression. All
methods can be used with covariate adjusted permutation. If the user doesn't
supply a covariate file, or disables covariate adjustment, the phenotypes
are shuffled uniformly.

## Supported Methods ##

Methods can be chosen using the -m, or --method option. The default method is 
VAAST.

	- BURDEN (Wu, Guan, Pankow 2016)
	- CALPHA (Neale et al. 2011)
	- CMC (Li, Leal 2008)
	- RVT1 (Morris, Zeggini 2010)
	- RVT2 (Morris, Zeggini 2010)
	- SKAT (Wu et al. 2011; Wu, Guan, Pankow 2016)
	- SKATO (Lee, Wu, Lin 2012; Wu, Guan, Pankow 2016)
	- VAAST (Yandell et al. 2011) -- Default
	- VT (Price et al. 2010)
	- WSS (Madsen, Browning 2009)

### SKAT / SKAT-O ###

SKAT and SKATO are implemented using the method described by Wu, Guan,
and Pankow (2016).  The variant based statistic provides significant
computational speedup over the individual based statistic
for very large datasets. Additionally, this allows SKAT to be used on very
large datasets, where otherwise the NxN covariance matrix will exhaust available
 memory. Even with the computation advantage that this
approach provides, it is still time consuming to permute SKAT-O.

Note that SKAT can return two different values. When used without permutation,
SKAT returns an analytic p-value. When used with permutation, the costly 
calculation of the p-value is skipped, and the test statistic is returned
instead.

## Compiling ##

Dependencies: 
- Tested with Armadillo >= 8.600
- C++ compiler supporting C++14
- C++ Boost Library > 1.66 (required for quadrature)

The dependencies can be installed via your package manager. Using Homebrew on
macOS:

```bash
brew install armadillo
brew install boost

```

If you're working in a cluster or otherwise managed environment,
ensure that Armadillo is compiled against both Lapack and BLAS, or an
alternative like Intel MKL. Without that, this software will fail when
calling Singular Value Decomposition and other algorithms not provided
by Armadillo itself.

Create a build directory and run cmake.

```bash
mkdir build
cd build

cmake -DCMAKE_BUILD_TYPE=Release ..
make
```

If cmake fails to detect armadillo, but you're sure it is available,
you may need to direct cmake to the library, e.g., when compiling on a
cluster, with packages in non-standard locations. In that case the
following should work:

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

The location for boost may need to be specified if it isn't installed in a
typical location.

```bash
cmake -DBOOST_ROOT=<path_to_boost> ..
```

You can combine the above as necessary. Earlier versions of the
Armadillo library may work, but haven't been tested. If you need to
change the compiler used from the one automatically detected to
another, perhaps newer compiler:

```bash
cmake -DCMAKE_CXX_COMPILER=<path_to_executable> ..
```

## Benchmarks ##

Testing on a mid-2015 MacBook Pro, 2.2Ghz Core i7, 16GB 1600 Mhz DDR3 ram.

| Test  | Samples | Genes | Permutations |  Time (sec) |
| -----:|:-------:|:-----:|:------------:|:-----------:|
| CMC   | 100000  | 1000  | 0            |  341.42     |
| RVT1  | 100000  | 1000  | 0            |  377.89     |
| RVT2  | 100000  | 1000  | 0            |  405.56     |
| VAAST | 100000  | 1000  | 10000        |  2086.53    |
| SKAT-O| 100000  | 1000  | 0            |  487.79     |
| WSS   | 100000  | 1000  | 10000        |  3432.29    |
| SKAT  | 100000  | 1000  | 10000        |  28293.36   |

### Gzipped Data ###

1049.7s with gzipped data, 1000 permutations of VAAST.
1063.4s with unzipped data, 1000 permutations of VAAST.

Working with gzipped data has no significant impact on runtime.

## File Formats ##

The file formats used are simple, plain text formats, which may be gzipped.

### Matrix ###

Matrix files may be zipped with gzip, or unzipped. They are provided
to the program with the "-i" flag. The matrix format used for the
genotype file includes a header and is is as follows:

	1) Chromosome
	2) Start position (bp)
	3) End position (bp)
	4) Type (SNV) 
	5) Reference allele
	6) Alternate allele
	7) Gene (e.g., BRCA1)
	8) Transcript (e.g., NM_700030) 
	9) Region (e.g., exonic, intronic)
	10) Function (synonymous, nonsynonymous)
	11) CASM score (function and conservation based weight)
	12+) Sample genotype (e.g., 0/1/2 - will be converted to minor allele count)

Note that the type field is used for
grouping in VAAST, and the program will fail if the
annotation is not present when grouping is used.

Variants are expected to be repeated for each transcript they appear
in. This does increase the size of the file, but simplifies
parsing. Because the format is simple, and QC is expected to be
finished before this program is run, the file remains reasonably
lean. E.g., with 20,000 simulated genes, and 100,000 samples, a matrix
file, uncompressed, is only 56GB. Compressed with gzip defaults, the same
matrix file is only 841MB. 

### Bed Mask File ###

The "-b" option allows the user to provide a ".bed" format file to mask
problematic variants. The .bed format for the mask file is as follows:

	1) Chromosome
	2) Start position of masked region
	3) End position of masked region
	
XQC (in the Cross Platform Association Toolkit (XPAT)) will provide a mask
file for those variants that fail QC if it is being used.

### Covariates ###

Covariates can be provided in a matrix format file. The file format is as
follows:

    1) Sample ID
    2) First covariate
    3) Second covariate
    4) ...

There is no limit on the number of covariates that can be provided but they are
assumed to be numeric.

### Weights ###

Essentially a bed format file, with the 4th column the variant annotation and 
5th column the weight, which will be converted to a log scale for VAAST. Weights
will be used by VAAST, and SKAT / SKAT-O if the linear kernel is used.

    1) Chromosome
    2) Start position
    3) End position
    4) Annotation (e.g. SNV, Deletion, Insertion, etc.)
    5) Weight

Weights provided via this option overwrite the CASM scores, which are used for
weighting by default. Passing the --no_weights option will remove all weights.
Using SKAT or SKAT-O with the weighted Linear kernel option will overwrite the
weights with those calculated from a beta distribution.

## Filter Whitelist ##

This is a list of variant types and functions to include for each method. Any
types or functions not in the list will be ignored / removed during parsing. The
whitelist format is simple. It is the method matched with the type / function.
Example:

- BURDEN:SNV
- BURDEN:insertion
- BURDEN:deletion
- BURDEN:SPDA
- BURDEN:complex_substitution
- BURDEN:synonymous
- BURDEN:nonsynonymous
- BURDEN:frameshift
- BURDEN:nonframeshift
- BURDEN:stopgain
- BURDEN:stoploss
- BURDEN:unknown

This is the default allowed set of variant types and functions for the BURDEN
test. Removing the SPDA line would result in all variants labelled as splicing
donor / acceptor being removed.