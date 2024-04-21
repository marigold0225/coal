cmake_minimum_required(VERSION 3.27)

find_package(benchmark REQUIRED)

set(BENCHMARK_SOURCES
        #        benchmark_test.cpp
                benchmark_lorentzMatrix.cpp)
#                benchmark_Jacobi.cpp)
#        benchmark_checkValidCombination.cpp)

add_executable(benchmark_test ${BENCHMARK_SOURCES})

target_link_libraries(benchmark_test benchmark::benchmark coal_lib)

target_include_directories(benchmark_test PRIVATE ../src/headers)

find_package(Threads REQUIRED)
target_link_libraries(benchmark_test Threads::Threads)
