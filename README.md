# Permute Associate

Permute Associate is a tool for gene-based testing of rare variants in sequencing
studies of any size. We provide methods for the analysis of data provided in simple,
plain text, array formats. Permutation and statistical computation is multithreaded.

This software has been carefully designed to reduce memory overhead as much as possible
through the use of shared read only objects, and elimination of data that is no longer
needed. Total concurrent memory usage increases with the number of threads. By default,
the software uses half of the available number of cpus, but this can be controlled.

All methods are available for the analysis of binary traits. Quantitative traits can be analyzed
using the BURDEN, RVT1, RVT2, SKAT, and SKATO methods currently.

## Approach

We have implemented an alternative to the efficient resampling approach
presented in Lee, Fuchsberger, Kim, and Scott (2016). In short, we use 
covariate adjusted permutation for the minor allele carriers in a gene. 
We permute the phenotype of all individuals according to their odds estimated from
logistic regression. This only applies to methods that do not use residuals from
a linear or glm fit. For other methods, phenotypes are uniformly shuffled.

SKAT and SKATO are implemented using the method described by Wu, Guan, and Pankow (2016).
The variant based statistic provides significant computational speedup over the individual
based statistic, especially for very large datasets. Even with the computation advantage that
this approach provides, it is still time consuming to permute SKAT-O.

## Supported Methods

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

## Compiling

Dependencies: 
- Tested with Armadillo >= 8.600
- C++ compiler supporting C++14
- C++ Boost Library > 1.66 (required for quadrature)

The dependencies can be installed via your package manager. Using Homebrew on macOS:

```bash
brew install armadillo
brew install boost

```

If you're working in a cluster or otherwise managed environment, ensure that Armadillo is compiled
against both Lapack and BLAS, or an alternative like Intel MKL. Without that, this software will fail
when calling Singular Value Decomposition and other algorithms not provided by Armadillo itself.

Create a build directory and run cmake.

```bash
mkdir pabuild
cd pabuild

cmake -DCMAKE_BUILD_TYPE=Release ..
make
```

If cmake fails to detect armadillo, but you're sure it is available, you may 
need to direct cmake to the library, e.g., when compiling on a cluster, with 
packages in non-standard locations. In that case the following should work:

```bash
mkdir pabuild && cd pabuild

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

You can combine the above as necessary. Earlier versions of the Armadillo
library may work, but haven't been tested. If you need to change the compiler used
from the one automatically detected to another, perhaps newer compiler:

```bash
cmake -DCMAKE_CXX_COMPILER=<path_to_executable> ..
```