cmake_minimum_required(VERSION 3.15)
project(CAPER VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(Boost_USE_MULTITHREADED TRUE)

option(ENABLE_NATIVE_OPTIMIZATIONS "Enable native architecture-specific optimizations" OFF)

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    message(STATUS "GNU C++ compiler detected")
endif()

find_package(Armadillo 8.6 REQUIRED)
find_package(Threads REQUIRED)
find_package(Boost 1.66.0 COMPONENTS program_options iostreams REQUIRED)
find_package(ZLIB REQUIRED)

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

if (UNIX)
    find_package(LAPACK REQUIRED)
    find_package(BLAS REQUIRED)
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


add_executable(caper caper/caper.cpp ${PROJECT_SUPPORT_FILES} utility/jointhreads.hpp)

if (FALSE)
    target_link_libraries(caper nvblas)
endif ()
target_include_directories(caper PRIVATE
        "${CMAKE_CURRENT_LIST_DIR}/third_party"
        "${CMAKE_CURRENT_LIST_DIR}/third_party/sse2neon")
target_link_libraries(caper PRIVATE
        Armadillo::armadillo
        Threads::Threads
        Boost::program_options
        Boost::iostreams
        LAPACK::LAPACK
        BLAS::BLAS
        ZLIB::ZLIB)

target_compile_options(caper PRIVATE
        "$<$<BOOL:${ENABLE_NATIVE_OPTIMIZATIONS}>:$<$<OR:$<CXX_COMPILER_ID:GNU>,$<CXX_COMPILER_ID:Clang>>:-march=native>>"
        "$<$<AND:$<BOOL:${ENABLE_NATIVE_OPTIMIZATIONS}>,$<CXX_COMPILER_ID:GNU>>:-fno-openmp>")

target_compile_definitions(caper PRIVATE MAXCOLORS=10000000)

add_executable(power power/power.cpp ${PROJECT_SUPPORT_FILES} utility/jointhreads.hpp)
target_include_directories(power PRIVATE
        "${CMAKE_CURRENT_LIST_DIR}/third_party"
        "${CMAKE_CURRENT_LIST_DIR}/third_party/sse2neon")
target_link_libraries(power PRIVATE
        Armadillo::armadillo
        Threads::Threads
        Boost::program_options
        Boost::iostreams
        LAPACK::LAPACK
        BLAS::BLAS)

target_compile_options(power PRIVATE
        "$<$<BOOL:${ENABLE_NATIVE_OPTIMIZATIONS}>:$<$<OR:$<CXX_COMPILER_ID:GNU>,$<CXX_COMPILER_ID:Clang>>:-march=native>>"
        "$<$<AND:$<BOOL:${ENABLE_NATIVE_OPTIMIZATIONS}>,$<CXX_COMPILER_ID:GNU>>:-fno-openmp>")

target_compile_definitions(power PRIVATE MAXCOLORS=10000000)

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
        ${PROJECT_SUPPORT_FILES})
target_include_directories(catch_test PRIVATE
        "${CMAKE_CURRENT_LIST_DIR}/third_party"
        "${CMAKE_CURRENT_LIST_DIR}/third_party/sse2neon")
target_link_libraries(catch_test PRIVATE
        Armadillo::armadillo
        Threads::Threads
        Boost::program_options
        Boost::iostreams
        LAPACK::LAPACK
        BLAS::BLAS)

target_compile_options(catch_test PRIVATE
        "$<$<BOOL:${ENABLE_NATIVE_OPTIMIZATIONS}>:$<$<OR:$<CXX_COMPILER_ID:GNU>,$<CXX_COMPILER_ID:Clang>>:-march=native>>"
        "$<$<AND:$<BOOL:${ENABLE_NATIVE_OPTIMIZATIONS}>,$<CXX_COMPILER_ID:GNU>>:-fno-openmp>")

target_compile_definitions(catch_test PRIVATE MAXCOLORS=10000000)
