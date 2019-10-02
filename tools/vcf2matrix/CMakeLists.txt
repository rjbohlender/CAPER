cmake_minimum_required(VERSION 2.8)
project(Permute_Associate C CXX)

set(CMAKE_CXX_STANDARD 14)

find_package(Boost 1.66.0 COMPONENTS program_options iostreams REQUIRED)

include_directories(${Boost_INCLUDE_DIRS})

add_executable(vcf2matrix
               vcf2matrix.cpp reference.cpp reference.hpp ../../utility/filesystem.cpp ../../utility/filesystem.hpp parser.cpp parser.hpp ../../utility/jointhreads.hpp)

target_link_libraries(vcf2matrix ${Boost_LIBRARIES})
