cmake_minimum_required(VERSION 3.10)
project(proName)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 不需要 set(CMAKE_PREFIX_PATH ...)，vcpkg 会通过 toolchain 自动处理

find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(ArrayFire REQUIRED)

message(STATUS "ArrayFire_INCLUDE_DIRS: ${ArrayFire_INCLUDE_DIRS}")
message(STATUS "ArrayFire_LIBRARIES: ${ArrayFire_LIBRARIES}")

add_executable(proName main.cpp)
target_include_directories(proName PRIVATE ${EIGEN3_INCLUDE_DIR})
target_include_directories(proName PRIVATE ${ArrayFire_INCLUDE_DIRS})

target_link_libraries(proName PRIVATE
        "C:/Users/29436/vcpkg/installed/x64-windows/lib/af.lib"
)


message(STATUS "ArrayFire_INCLUDE_DIRS: ${ArrayFire_INCLUDE_DIRS}")
message(STATUS "ArrayFire_LIBRARIES: ${ArrayFire_LIBRARIES}")

add_executable(testEigen textforEigen.cpp)
target_include_directories(testEigen PRIVATE ${EIGEN3_INCLUDE_DIR})
