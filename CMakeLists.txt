cmake_minimum_required(VERSION 2.8)
project(Permute_Associate)

set(CMAKE_CXX_STANDARD 14)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

find_package(Armadillo REQUIRED)
find_package(Threads REQUIRED)
find_package(Boost 1.66.0 COMPONENTS program_options REQUIRED)

if (APPLE)
    # -fsanitize=address -fsanitize=undefined -Wall -pedantic
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -pg")
endif ()

if (UNIX)
    find_package(LAPACK REQUIRED)
    find_package(BLAS REQUIRED)
endif()

if (CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_CXX_FLAGS_RELEASE  "${CMAKE_CXX_FLAGS_RELEASE} -fno-openmp")
else()
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=native")
endif ()

message(${CMAKE_CXX_FLAGS_DEBUG})
message(${CMAKE_CXX_FLAGS_RELEASE})

set(PROJECT_SUPPORT_FILES
    data/covariates.hpp
    utility/split.hpp
    statistics/methods.cpp
    statistics/methods.hpp
    data/gene.hpp
    data/taskargs.hpp
    data/taskqueue.hpp
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
    data/taskargs.cpp
    data/taskqueue.cpp
    statistics/skat.cpp
    statistics/skat.hpp
    statistics/skat_adjust.cpp
    statistics/skat_adjust.hpp
    data/covariates.cpp
    utility/main_support.cpp
    utility/main_support.hpp
    statistics/vaast.cpp
    statistics/vaast.hpp
    utility/jobdispatcher.cpp
    utility/jobdispatcher.hpp
    utility/reporter.cpp
    utility/reporter.hpp
    data/taskparams.hpp
    link/gaussian.cpp
    link/gaussian.hpp
    link/family.cpp
    link/family.hpp
    link/binomial.cpp
    link/binomial.hpp
    statistics/glm.hpp
    statistics/bayesianglm.hpp statistics/vt.cpp statistics/vt.hpp)

include_directories("${CMAKE_CURRENT_LIST_DIR}/third_party")
include_directories(${ARMADILLO_INCLUDE_DIRS})
message(${ARMADILLO_INCLUDE_DIRS})
include_directories(${Boost_INCLUDE_DIRS})
include_directories(${LAPACK_INCLUDE_DIRS})
include_directories(${BLAS_INCLUDE_DIRS})

add_executable(Permute_Associate
               main.cpp
               ${PROJECT_SUPPORT_FILES} data/calccasm.cpp data/calccasm.hpp statistics/skatr.cpp statistics/skatr.hpp)

target_link_libraries(Permute_Associate ${ARMADILLO_LIBRARIES})
target_link_libraries(Permute_Associate ${CMAKE_THREAD_LIBS_INIT})
target_link_libraries(Permute_Associate ${Boost_LIBRARIES})
target_link_libraries(Permute_Associate ${LAPACK_LIBRARIES})
target_link_libraries(Permute_Associate ${BLAS_LIBRARIES})

target_compile_definitions(Permute_Associate PRIVATE MAXCOLORS=10000000)

add_executable(catch_test
               catch_test/test.cpp
               ${PROJECT_SUPPORT_FILES})

target_include_directories(catch_test PRIVATE "${CMAKE_CURRENT_LIST_DIR}/third_party")
target_include_directories(catch_test PRIVATE ${ARMADILLO_INCLUDE_DIR})

target_link_libraries(catch_test ${ARMADILLO_LIBRARIES})
target_link_libraries(catch_test ${CMAKE_THREAD_LIBS_INIT})
