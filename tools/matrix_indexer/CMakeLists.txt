cmake_minimum_required(VERSION 3.15)
project(matrix_indexer VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)

find_package(Boost 1.66.0 COMPONENTS program_options iostreams REQUIRED)

include_directories(${Boost_INCLUDE_DIRS})

add_executable(matrix_indexer matrix_indexer.cpp)

target_link_libraries(matrix_indexer ${Boost_LIBRARIES})
