// Copyright Contributors to the OpenVDB Project
// SPDX-License-Identifier: MPL-2.0

#include "FlipExample.h"

// TO BUILD:
// mkdir build
// cd build
// cmake -DOPENVDB_BUILD_EXAMPLES=ON ../
// make -j 8
int
main(int argc, char *argv[])
{
    openvdb::initialize();

    example::FlipSolver flipSim(0.1f /* voxel size */);
    flipSim.render();
}