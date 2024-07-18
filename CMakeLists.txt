cmake_minimum_required(VERSION 3.27)

project(coal)

set(CMAKE_CXX_STANDARD 20)
# #

if (NOT CMAKE_BUILD_TYPE)
        set(CMAKE_BUILD_TYPE Release)
endif ()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -fopenmp -fno-math-errno -march=native -ftree-vectorize -mavx -mfma")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -fopenmp -fno-math-errno -march=native -ftree-vectorize -mavx -mfma")

set(LIB_SOURCE_FILES
        src/headers/Centrality.h
        src/cpp/Centrality.cpp
        src/headers/method.h
        src/cpp/method.cpp
        src/headers/particle.h
        src/cpp/particle.cpp
        src/headers/run.h
        src/cpp/run.cpp
        src/headers/flow.h
        src/cpp/flow.cpp
        src/headers/random.h
        src/headers/fileloader.h
        src/cpp/fileloader.cpp
        src/headers/chrono.h
        src/headers/output.h
        src/headers/config.h
        src/cpp/config.cpp
        src/headers/shortcut.h
        src/headers/ClusterCalculator.h
        src/headers/ThreadPool.h
        src/cpp/ThreadPool.cpp
        src/cpp/ClusterCalculator.cpp
        src/headers/logger.h
        src/headers/MatrixMethod.h
        src/headers/ctpl.h
        src/headers/xmath.h
        src/cpp/xmath.cpp)


#find_package(TBB REQUIRED)
# lib
add_library(coal_lib STATIC ${LIB_SOURCE_FILES})

# exec
add_executable(${PROJECT_NAME} src/cpp/main.cpp)

target_link_libraries(${PROJECT_NAME} coal_lib)
#target_link_libraries(coal_lib TBB::tbb)

# GSL
#find_package(GSL REQUIRED)
#if (NOT GSL_FOUND)
#    message(FATAL_ERROR "GSL not found. Please specify the GSL installation path using -DGSL_DIR=<path> or -DCMAKE_PREFIX_PATH=<path>.")
#endif ()
#include_directories(${GSL_INCLUDE_DIRS})
#target_link_libraries(coal_lib ${GSL_LIBRARIES})

# spdlog
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/3rdparty/spdlog")
        #    add_subdirectory(3rdparty/spdlog)
        include_directories(3rdparty/spdlog/include)
else()
        message(FATAL_ERROR "spdlog directory not found in 3rdparty. Please download spdlog and extract it to ${CMAKE_CURRENT_SOURCE_DIR}/3rdparty/spdlog")
endif()

# YAML
if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/3rdparty/yaml-cpp")
        add_subdirectory(3rdparty/yaml-cpp)
        target_link_libraries(coal_lib yaml-cpp)
else ()
        message(FATAL_ERROR "yaml-cpp directory not found in 3rdparty. Please download yaml-cpp from https://github.com/jbeder/yaml-cpp and extract it to ${CMAKE_CURRENT_SOURCE_DIR}/3rdparty/yaml-cpp")
endif ()

# Eigen
if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/3rdparty/eigen-3.4.0")
        include_directories(3rdparty/eigen-3.4.0)
else ()
        message(FATAL_ERROR "Eigen directory not found in 3rdparty. Please download Eigen from https://gitlab.com/libeigen/eigen/-/releases/3.4.0 and extract it to ${CMAKE_CURRENT_SOURCE_DIR}/3rdparty/eigen-3.4.0")
endif ()


# MKL
#include_directories(/home/mafu/intel/oneapi/mkl/latest/include)
#link_directories(/home/mafu/intel/oneapi/mkl/2024.1/lib /home/mafu/intel/oneapi/compiler/latest/lib)
#target_link_libraries(coal_lib
#        /home/mafu/intel/oneapi/mkl/2024.1/lib/libmkl_rt.so
#)


# Google Test
option(PACKAGE_TESTS "Build the tests" OFF)
if (PACKAGE_TESTS AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/googletest")
        add_subdirectory(googletest)
        include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})
        add_subdirectory(tests)
else ()
        message("Google Test directory not found, skipping tests")
endif ()

# benchmark
OPTION(BUILD_BENCHMARK "Builds the benchmark subproject" OFF)
if (BUILD_BENCHMARK)
        find_package(benchmark)
        if (benchmark_FOUND)
                message("BENCHMARK FOUND")
                add_subdirectory(benchmark)
        else ()
                message("Benchmark library not found, skipping benchmarks")
        endif ()
endif ()





