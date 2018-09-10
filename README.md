# Permute Associate

## Approach

We have implemented an alternative to the efficient resampling approach
presented in Lee, Fuchsberger, Kim, and Scott (2016). 

## Supported Methods

Methods can be chosen using the -m, or --method option. The default method is 
VAAST.

- CALPHA
- CMC
- SKAT
- VAAST
- VT
- WSS

## Compiling

Dependencies: 
- Tested with Armadillo 8.600
- C++ compiler supporting C++14
- C++ Boost Library

The dependencies can be installed via your package manager. Using Homebrew on macOS:

```bash
brew install armadillo
brew install boost

```

Create a build directory and run cmake.

```bash
mkdir pabuild
cd pabuild

cmake -DCMAKE_BUILD_TYPE=Release ..
make
```

If cmake fails to detect armadillo, but you're sure it is available, you may 
need to direct cmake to the library, e.g., when compiling on a cluster. In that
case the following should work:

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

You can combine the above as necessary. Earlier versions of the Armadillo
library may work, but haven't been tested. If you need to change the compiler used
from the one automatically detected to another, perhaps newer compiler:

```bash
cmake -DCMAKE_CXX_COMPILER=<path_to_executable> ..
```
