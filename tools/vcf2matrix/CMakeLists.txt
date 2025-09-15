cmake_minimum_required(VERSION 3.15)
project(vcf2matrix VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_executable(vcf2matrix
    vcf2matrix.cpp
    parser.cpp
    reference.cpp)

target_include_directories(vcf2matrix PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(vcf2matrix PRIVATE
    caper_utils
    Boost::program_options
    Boost::iostreams)

target_compile_options(vcf2matrix PRIVATE
    "$<$<BOOL:${ENABLE_NATIVE_OPTIMIZATIONS}>:$<$<OR:$<CXX_COMPILER_ID:GNU>,$<CXX_COMPILER_ID:Clang>>:-march=native>>"
    "$<$<AND:$<BOOL:${ENABLE_NATIVE_OPTIMIZATIONS}>,$<CXX_COMPILER_ID:GNU>>:-fno-openmp>")

install(TARGETS vcf2matrix RUNTIME DESTINATION bin)

