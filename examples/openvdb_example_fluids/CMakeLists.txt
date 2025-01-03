# Copyright Contributors to the OpenVDB Project
# SPDX-License-Identifier: Apache-2.0
#
#[=======================================================================[
    CMake Configuration for OpenVDB Example Fluids
#]=======================================================================]

cmake_minimum_required(VERSION 3.18)
project(OpenVDBExampleFluids LANGUAGES CXX)

include(GNUInstallDirs)

set(SOURCE_FILES main.cc)
add_executable(openvdb_example_fluids ${SOURCE_FILES})
target_link_libraries(openvdb_example_fluids ${OPENVDB_BINARIES_DEPENDENT_LIBS})

install(TARGETS openvdb_example_fluids RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})