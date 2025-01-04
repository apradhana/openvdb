// Copyright Contributors to the OpenVDB Project
// SPDX-License-Identifier: Apache-2.0

#include "FlipExample.h"
#include "SmokeExample.h"

// TO BUILD:
// mkdir build
// cd build
// cmake -DOPENVDB_BUILD_EXAMPLES=ON ../
// make -j 8
int
main(int argc, char *argv[])
{
    openvdb::initialize();

#if 0
    example::FlipSolver flipSim(0.1f /* voxel size */);
    flipSim.render();
#endif
#if 1
    example::SmokeSolver smokeSim(0.1f /* voxel size */);
    smokeSim.render();
#endif
}