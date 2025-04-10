cmake_minimum_required(VERSION 2.8)
project(Permute_Associate C CXX)

set(CMAKE_CXX_STANDARD 14)

find_package(Boost 1.66.0 COMPONENTS program_options iostreams REQUIRED)

include_directories(${Boost_INCLUDE_DIRS})

add_executable(matrix_indexer matrix_indexer.cpp)

target_link_libraries(matrix_indexer ${Boost_LIBRARIES})
