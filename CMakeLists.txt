cmake_minimum_required(VERSION 3.15)
project(CAPER VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include(FetchContent)

set(Boost_USE_MULTITHREADED TRUE)

option(ENABLE_NATIVE_OPTIMIZATIONS "Enable native architecture-specific optimizations" OFF)

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    message(STATUS "GNU C++ compiler detected")
endif()

find_package(Armadillo 12.6 QUIET)
set(USING_SYSTEM_ARMADILLO ${Armadillo_FOUND})

if (NOT USING_SYSTEM_ARMADILLO)
    message(STATUS "System Armadillo >= 12.6 not found, fetching 12.6.4")
    include(FetchContent)
    FetchContent_Declare(armadillo
        URL https://downloads.sourceforge.net/project/arma/armadillo-12.6.4.tar.xz
        URL_HASH SHA256=eb7f243ffc32f18324bc7fa978d0358637e7357ca7836bec55b4eb56e9749380
        DOWNLOAD_EXTRACT_TIMESTAMP TRUE)
    FetchContent_MakeAvailable(armadillo)

    if (NOT TARGET Armadillo::armadillo)
        add_library(Armadillo::armadillo INTERFACE IMPORTED)
    endif ()

    set(ARMADILLO_INCLUDE_DIRS "${armadillo_SOURCE_DIR}/include")
    set(ARMADILLO_LIBRARIES "")
    set_target_properties(Armadillo::armadillo PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${ARMADILLO_INCLUDE_DIRS}")
    target_compile_definitions(Armadillo::armadillo INTERFACE ARMA_DONT_USE_WRAPPER ARMA_DONT_USE_LAPACK ARMA_DONT_USE_BLAS)
else ()
    message(STATUS "Using system Armadillo ${Armadillo_VERSION}")
endif ()

find_package(Threads REQUIRED)
find_package(Boost 1.66.0 COMPONENTS program_options iostreams QUIET)
set(USING_SYSTEM_BOOST ${Boost_FOUND})

if (NOT Boost_FOUND)
    message(STATUS "System Boost >= 1.66 not found, fetching 1.82.0")
    FetchContent_Declare(boost
        URL https://downloads.sourceforge.net/project/boost/boost/1.82.0/boost_1_82_0.tar.gz
        URL_HASH SHA256=66a469b6e608a51f8347236f4912e27dc5c60c60d7d53ae9bfe4683316c6f04c
        DOWNLOAD_EXTRACT_TIMESTAMP TRUE)
    FetchContent_MakeAvailable(boost)

    FetchContent_Declare(zstd
        URL https://github.com/facebook/zstd/archive/refs/tags/v1.5.5.tar.gz
        URL_HASH SHA256=98e9c3d949d1b924e28e01eccb7deed865eefebf25c2f21c702e5cd5b63b85e1
        DOWNLOAD_EXTRACT_TIMESTAMP TRUE
        SOURCE_SUBDIR build/cmake)
    FetchContent_MakeAvailable(zstd)

    set(Boost_FOUND TRUE)
    set(Boost_VERSION 108200)
    set(Boost_INCLUDE_DIRS "${boost_SOURCE_DIR}")

    file(GLOB BOOST_PROGRAM_OPTIONS_SOURCES
        "${boost_SOURCE_DIR}/libs/program_options/src/*.cpp")

    add_library(boost_program_options STATIC ${BOOST_PROGRAM_OPTIONS_SOURCES})
    target_include_directories(boost_program_options PUBLIC "${Boost_INCLUDE_DIRS}")
    target_compile_definitions(boost_program_options PUBLIC BOOST_ALL_NO_LIB)

    add_library(boost_iostreams STATIC
        "${boost_SOURCE_DIR}/libs/iostreams/src/file_descriptor.cpp"
        "${boost_SOURCE_DIR}/libs/iostreams/src/gzip.cpp"
        "${boost_SOURCE_DIR}/libs/iostreams/src/mapped_file.cpp"
        "${boost_SOURCE_DIR}/libs/iostreams/src/zlib.cpp"
        "${boost_SOURCE_DIR}/libs/iostreams/src/zstd.cpp")
    target_include_directories(boost_iostreams PUBLIC "${Boost_INCLUDE_DIRS}")
    target_compile_definitions(boost_iostreams PUBLIC BOOST_ALL_NO_LIB BOOST_IOSTREAMS_NO_BZIP2 BOOST_IOSTREAMS_NO_LZMA)
    target_link_libraries(boost_iostreams PUBLIC ZLIB::ZLIB libzstd_static)

    add_library(Boost::program_options ALIAS boost_program_options)
    add_library(Boost::iostreams ALIAS boost_iostreams)
endif ()
find_package(ZLIB REQUIRED)

add_library(caper_utils STATIC utility/filesystem.cpp)
target_include_directories(caper_utils PUBLIC
        "${CMAKE_CURRENT_LIST_DIR}")
target_link_libraries(caper_utils PUBLIC Boost::iostreams)

if(NOT TARGET Armadillo::armadillo)
    add_library(Armadillo::armadillo INTERFACE IMPORTED)
    set_target_properties(Armadillo::armadillo PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${ARMADILLO_INCLUDE_DIRS}"
        INTERFACE_LINK_LIBRARIES "${ARMADILLO_LIBRARIES}")
endif()

if (APPLE)
    # -fsanitize=address -fsanitize=undefined -Wall -pedantic
    # set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O1 -fsanitize=address -fsanitize-blacklist=/Users/rjbohlender/CLionProjects/CAPER/Asan.supp -Wall -pedantic")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -pedantic")
endif ()

if (USING_SYSTEM_ARMADILLO)
    if (UNIX)
        find_package(LAPACK REQUIRED)
        find_package(BLAS REQUIRED)
    endif ()
endif ()

add_subdirectory(tools/vcf2matrix)
add_subdirectory(tools/matrix_indexer)

set(PROJECT_SUPPORT_FILES
    data/covariates.hpp
    utility/split.hpp
    statistics/methods.cpp
    statistics/methods.hpp
    data/gene.hpp
    caper/capertask.hpp
    third_party/stocc/stocc.h
    third_party/stocc/stoc1.cpp
    third_party/stocc/stoc2.cpp
    third_party/stocc/stoc3.cpp
    third_party/stocc/fnchyppr.cpp
    third_party/stocc/wnchyppr.cpp
    third_party/stocc/erfres.cpp
    third_party/stocc/randomc.h
    third_party/stocc/mersenne.cpp
    third_party/stocc/userintf.cpp
    third_party/stocc/sfmt.h
    third_party/stocc/sfmt.cpp
    statistics/logistic_regression.hpp
    statistics/logistic_regression.cpp
    data/permutation.cpp
    data/permutation.hpp
    data/bed.hpp
    data/weight.hpp
    data/result.hpp
    utility/filesystem.cpp
    utility/filesystem.hpp
    data/gene.cpp
    data/bed.cpp
    data/result.cpp
    data/weight.cpp
    third_party/QFC/qfc2.cpp
    third_party/QFC/qfc2.hpp
    caper/capertask.cpp
    data/covariates.cpp
    statistics/vaast.cpp
    statistics/vaast.hpp
    utility/jobdispatcher.hpp
    utility/reporter.cpp
    utility/reporter.hpp
    utility/taskparams.hpp
    link/gaussian.cpp
    link/gaussian.hpp
    link/family.cpp
    link/family.hpp
    link/binomial.cpp
    link/binomial.hpp
    statistics/glm.hpp
    statistics/bayesianglm.hpp
    statistics/vt.cpp
    statistics/vt.hpp
    statistics/skat.cpp
    statistics/skat.hpp
    power/powerop.cpp
    power/powerop.hpp
    utility/taskqueue.hpp
    caper/caperop.cpp
    caper/caperop.hpp
    power/powertask.cpp
    power/powertask.hpp
    caese/caeseop.cpp
    caese/caeseop.hpp
    caese/caesetask.cpp
    caese/caesetask.hpp
    statistics/fishertest.cpp
    statistics/fishertest.hpp
    utility/jointhreads.hpp
    utility/indexsort.hpp
    utility/math.hpp
    utility/filevalidator.cpp
    utility/filevalidator.hpp
    data/filter.cpp
    data/filter.hpp
    data/matrix_indices.hpp)


add_library(caper_core STATIC ${PROJECT_SUPPORT_FILES})
target_include_directories(caper_core PUBLIC
        "${CMAKE_CURRENT_LIST_DIR}/third_party"
        "${CMAKE_CURRENT_LIST_DIR}/third_party/sse2neon")
target_link_libraries(caper_core PUBLIC
        Armadillo::armadillo
        Threads::Threads
        Boost::program_options
        Boost::iostreams
        ZLIB::ZLIB)

if (USING_SYSTEM_ARMADILLO)
    target_link_libraries(caper_core PUBLIC LAPACK::LAPACK BLAS::BLAS)
endif ()

target_compile_options(caper_core PUBLIC
        "$<$<BOOL:${ENABLE_NATIVE_OPTIMIZATIONS}>:$<$<OR:$<CXX_COMPILER_ID:GNU>,$<CXX_COMPILER_ID:Clang>>:-march=native>>"
        "$<$<AND:$<BOOL:${ENABLE_NATIVE_OPTIMIZATIONS}>,$<CXX_COMPILER_ID:GNU>>:-fno-openmp>")

target_compile_definitions(caper_core PUBLIC MAXCOLORS=10000000)


add_executable(caper caper/caper.cpp)

if (FALSE)
    target_link_libraries(caper nvblas)
endif ()
target_link_libraries(caper PRIVATE caper_core)

add_executable(power power/power.cpp)
target_link_libraries(power PRIVATE caper_core)

# add_executable(caese
#         caese/caese.cpp
#         ${PROJECT_SUPPORT_FILES} caese/caese.cpp utility/jointhreads.hpp)
#
# target_link_libraries(caese ${ARMADILLO_LIBRARIES})
# target_link_libraries(caese ${CMAKE_THREAD_LIBS_INIT})
# target_link_libraries(caese ${Boost_LIBRARIES})
# target_link_libraries(caese ${LAPACK_LIBRARIES})
# target_link_libraries(caese ${BLAS_LIBRARIES})
#
# target_compile_definitions(caese PRIVATE MAXCOLORS=10000000)

add_executable(catch_test
        catch_test/test.cpp
        catch_test/jobdispatcher_test.cpp
        catch_test/permutation_split_test.cpp
        catch_test/powerop_sample_test.cpp
        catch_test/permutation_test.cpp)
target_link_libraries(catch_test PRIVATE caper_core)
