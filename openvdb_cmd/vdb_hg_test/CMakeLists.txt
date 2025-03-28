# Copyright Contributors to the OpenVDB Project
# SPDX-License-Identifier: Apache-2.0
#
#[=======================================================================[

  CMake Configuration for VDB Half Grid Test

#]=======================================================================]

cmake_minimum_required(VERSION 3.20)
project(VDBHGTest LANGUAGES CXX)

include(GNUInstallDirs)

if(USE_IMATH_HALF)
  find_package(Imath ${MINIMUM_IMATH_VERSION} CONFIG REQUIRED)
endif()

set(SOURCE_FILES main.cc)
add_executable(vdb_hg_test ${SOURCE_FILES})

# Set deps. Note that the order here is important. If we're building against
# Houdini 17.5 we must include OpenEXR and IlmBase deps first to ensure the
# users chosen namespaced headers are correctly prioritized. Otherwise other
# include paths from shared installs (including houdini) may pull in the wrong
# headers

target_link_libraries(vdb_hg_test
  # For Imath/OpenEXR v3.X
  $<TARGET_NAME_IF_EXISTS:Imath::Imath>
  $<TARGET_NAME_IF_EXISTS:OpenEXR::IlmThread>
  $<TARGET_NAME_IF_EXISTS:OpenEXR::Iex>
  ${OPENVDB_BINARIES_DEPENDENT_LIBS}
)

install(TARGETS vdb_hg_test RUNTIME DESTINATION ${OPENVDB_INSTALL_BINDIR})
