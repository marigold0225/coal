cmake_minimum_required(VERSION 3.27)

project(coal)

set(CMAKE_CXX_STANDARD 20)

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Ofast")
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Ofast")
##
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native -ftree-vectorize -mavx -mfma")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -march=native -ftree-vectorize -mavx -mfma")

set(LIB_SOURCE_FILES
        src/headers/smash.h
        src/cpp/smash.cpp
        src/headers/method.h
        src/cpp/method.cpp
        src/headers/clusterMix.h
        src/cpp/clusterMix.cpp
        src/headers/particle.h
        src/cpp/particle.cpp
        src/headers/run.h
        src/cpp/run.cpp)


#find_package(TBB REQUIRED)


# lib
add_library(coal_lib STATIC ${LIB_SOURCE_FILES})

# exec
add_executable(${PROJECT_NAME} src/cpp/main.cpp)
#target_link_libraries(coal_lib TBB::tbb)
target_link_libraries(coal_lib)
target_link_libraries(${PROJECT_NAME} coal_lib)

# 3rd
add_subdirectory(3rdparty/yaml-cpp)
target_link_libraries(coal_lib yaml-cpp)

# openMP 配置
FIND_PACKAGE(OpenMP REQUIRED)
if (OPENMP_FOUND)
    message("OPENMP FOUND")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif ()

# Eigen
include_directories(3rdparty/eigen-3.4.0)

#Google Test
add_subdirectory(googletest)
include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

# test
#add_subdirectory(tests)





