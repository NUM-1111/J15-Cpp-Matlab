cmake_minimum_required(VERSION 3.30)
project(proName)

set(CMAKE_CXX_STANDARD 14)

set(CMAKE_PREFIX_PATH "C:/Users/29436/vcpkg/installed/x64-windows")


find_package(Eigen3 CONFIG REQUIRED)

add_executable(${PROJECT_NAME} main.cpp)

target_link_libraries(${PROJECT_NAME} PRIVATE Eigen3::Eigen)
