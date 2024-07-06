cmake_minimum_required(VERSION 3.27)
set(TEST_FILES
#        test_cluster.cpp
        test_method.cpp)

add_executable(runUnitTests ${TEST_FILES})

target_link_libraries(runUnitTests coal_lib gtest gtest_main)

target_include_directories(runUnitTests PRIVATE ../src/headers)