cmake_minimum_required(VERSION 3.10)
project(proName)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 指定 vcpkg 安装路径
set(CMAKE_PREFIX_PATH "C:/Users/29436/vcpkg/installed/x64-windows")

# 查找 Eigen
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

# 目标 1: 编译主文件 main.cpp
add_executable(proName main.cpp)
target_include_directories(proName PRIVATE ${EIGEN3_INCLUDE_DIR})

# 目标 2: 编译 textforEigen.cpp 作为单独的可执行文件
add_executable(testEigen textforEigen.cpp)
target_include_directories(testEigen PRIVATE ${EIGEN3_INCLUDE_DIR})
