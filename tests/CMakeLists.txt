# 设置最小CMake版本
cmake_minimum_required(VERSION 3.27)
# 定义测试文件
set(TEST_FILES
        test_cluster.cpp
        test_method.cpp)
# test_method.cpp)  # 添加更多的测试文件

# 添加测试可执行文件
add_executable(runUnitTests ${TEST_FILES})

# 链接到主项目库和Google Test
target_link_libraries(runUnitTests coal_lib gtest gtest_main)

# 包含主项目的头文件路径
target_include_directories(runUnitTests PRIVATE ../src/headers)