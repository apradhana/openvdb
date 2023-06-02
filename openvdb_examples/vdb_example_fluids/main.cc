// Copyright Contributors to the OpenVDB Project
// SPDX-License-Identifier: MPL-2.0

#include <algorithm>
#include <iostream>
#include <vector>
#include <openvdb/openvdb.h>
#include <openvdb/util/logging.h>


namespace openvdb_example_fluids {
class FluidSolver {
    private:
    virtual void initialize() = 0;

    protected:
    FluidSolver(float voxelSize)
        : mVoxelSize(voxelSize)
    {
        std::cout << "Fluid solver constructor" << std::endl;
    }

    ~FluidSolver()
    { std::cout << "FluidSolver destructor" << std::endl; }

    virtual void advanceOneSubstep(float const dt) = 0;

    virtual void writeState() = 0;

    protected:
        float mVoxelSize = 0.1f;
        // std::vector<openvdb::FloatGrid::Ptr> mColliders;
};


class FlipSolver : public FluidSolver {
    void initialize() {}

    public:

    FlipSolver(float const voxelSize)
        : FluidSolver{voxelSize}
    {
        std::cout << "FlipSolver constructor" << std::endl;
    }

    ~FlipSolver()
    { std::cout << "FlipSolver destructor" << std::endl; }

    void advanceOneSubstep(float const dt) override {
        std::cout << "FlipSolver::advanceOneSubstep dt = " << dt << std::endl;
    }

    void writeState() override {
        std::cout << "FlipSolver::writeState" << std::endl;
    }

    private:

};


class SmokeSolver : public FluidSolver {
    void initialize() {}

    public:

    SmokeSolver(float const voxelSize)
        : FluidSolver{voxelSize}
    {
        std::cout << "Smoke solver constructor" << std::endl;
    }

    ~SmokeSolver()
    { std::cout << "Smoke Solver destructor" << std::endl; }

    void advanceOneSubstep(float const dt) override {
        std::cout << "SmokeSolver::advanceOneSubstep dt = " << dt << std::endl;
    }

    void writeState() override {
        std::cout << "SmokeSolver::writeState" << std::endl;
    }

    private:
};

} // namespace openvdb_example_fluids

int
main(int argc, char *argv[])
{
    using namespace openvdb_example_fluids;

    openvdb::initialize();

    float const dt = 1.f / 24.f;
    int const numFrames = 30;
    float const voxelSize = 0.02;

    {
        FlipSolver solver(voxelSize);

        for (int i = 0; i < numFrames; i++) {
            solver.advanceOneSubstep(dt);
            solver.writeState();
        }
    }

    {
        SmokeSolver solver(voxelSize);

        for (int i = 0; i < numFrames; i++) {
            solver.advanceOneSubstep(dt);
            solver.writeState();
        }
    }
}
