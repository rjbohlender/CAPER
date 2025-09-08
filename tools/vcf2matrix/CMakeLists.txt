cmake_minimum_required(VERSION 3.15)
project(vcf2matrix VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)

find_package(Boost 1.66.0 COMPONENTS program_options iostreams REQUIRED)

include_directories(${Boost_INCLUDE_DIRS})

add_executable(vcf2matrix
               vcf2matrix.cpp reference.cpp reference.hpp ../../utility/filesystem.cpp ../../utility/filesystem.hpp parser.cpp parser.hpp ../../utility/jointhreads.hpp)

target_link_libraries(vcf2matrix ${Boost_LIBRARIES})

target_compile_options(vcf2matrix PRIVATE
        "$<$<BOOL:${ENABLE_NATIVE_OPTIMIZATIONS}>:$<$<OR:$<CXX_COMPILER_ID:GNU>,$<CXX_COMPILER_ID:Clang>>:-march=native>>"
        "$<$<AND:$<BOOL:${ENABLE_NATIVE_OPTIMIZATIONS}>,$<CXX_COMPILER_ID:GNU>>:-fno-openmp>")
