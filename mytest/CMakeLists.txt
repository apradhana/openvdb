cmake_minimum_required(VERSION 3.27)
project(test LANGUAGES CXX CUDA)

include_directories(
    /usr/local/cuda/include
)

set(CMAKE_CXX_STANDARD 17)

list(APPEND CMAKE_MODULE_PATH "/usr/local/lib/cmake/OpenVDB")
find_package(OpenVDB "11.0.1" REQUIRED)
add_executable(nanovdb_test main.cpp)

target_link_libraries(nanovdb_test OpenVDB::openvdb)
