cmake_minimum_required(VERSION 3.15)
project(matrix_indexer VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_executable(matrix_indexer matrix_indexer.cpp)

target_include_directories(matrix_indexer PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(matrix_indexer PRIVATE
    Boost::program_options
    Boost::iostreams)

target_compile_options(matrix_indexer PRIVATE
    "$<$<BOOL:${ENABLE_NATIVE_OPTIMIZATIONS}>:$<$<OR:$<CXX_COMPILER_ID:GNU>,$<CXX_COMPILER_ID:Clang>>:-march=native>>"
    "$<$<AND:$<BOOL:${ENABLE_NATIVE_OPTIMIZATIONS}>,$<CXX_COMPILER_ID:GNU>>:-fno-openmp>")

install(TARGETS matrix_indexer RUNTIME DESTINATION bin)

