// Copyright Contributors to the OpenVDB Project
// SPDX-License-Identifier: MPL-2.0

#include <algorithm>
#include <iostream>
#include <vector>
#include <openvdb/openvdb.h>
#include <openvdb/util/logging.h>


namespace openvdb_example_fluids {
class FluidSolver {
    void initialize() {}
    public:
    void advanceOneSubstep(float const dt) {
        std::cout << "FluidSolver::advanceOneSubstep dt = " << dt << std::endl;
    }

    void writeState() {
        std::cout << "FluidSolver::writeState" << std::endl;
    }
};


class FlipSolver {
    void initialize() {}
    public:
    void advanceOneSubstep(float const dt) {
        std::cout << "FlipSolver::advanceOneSubstep dt = " << dt << std::endl;
    }

    void writeState() {
        std::cout << "FlipSolver::writeState" << std::endl;
    }

};


class SmokeSolver {
    void initialize() {}
    public:
    void advanceOneSubstep(float const dt) {
        std::cout << "SmokeSolver::advanceOneSubstep dt = " << dt << std::endl;
    }

    void writeState() {
        std::cout << "SmokeSolver::writeState" << std::endl;
    }
};

} // namespace openvdb_example_fluids

int
main(int argc, char *argv[])
{
    using namespace openvdb_example_fluids;
    float const dt = 1.f / 24.f;
    int const numFrames = 30;
    {
        FlipSolver solver;
    
        for (int i = 0; i < numFrames; i++) {
            solver.advanceOneSubstep(dt);
            solver.writeState();
        }
    }

    { 
        SmokeSolver solver;
    
        for (int i = 0; i < numFrames; i++) {
            solver.advanceOneSubstep(dt);
            solver.writeState();
        }
    }


    std::cout << "Hello world!" << std::endl;
}