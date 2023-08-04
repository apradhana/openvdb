// Copyright Contributors to the OpenVDB Project
// SPDX-License-Identifier: MPL-2.0

# pragma once

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <math.h>
#include <string>

#include <openvdb/openvdb.h>
#include <openvdb/tools/Composite.h> // for tools::compMax
#include <openvdb/tools/GridOperators.h> // for divergence and gradient
#include <openvdb/tools/LevelSetSphere.h> // for createLevelSetSphere
#include <openvdb/tools/MeshToVolume.h> // for createLevelSetBox
#include <openvdb/tools/Morphology.h> // for erodeActiveValues
#include <openvdb/tools/PoissonSolver.h> // for poisson solve
#include <openvdb/tools/VolumeAdvect.h> // for tools::VolumeAdvection
#include <openvdb/tree/LeafManager.h> // for LeafManager

using namespace openvdb;

namespace example {

class SmokeSolver {

public:

    SmokeSolver(float const voxelSize);

    void render();

    void foobar();

private:

    void initialize();

    void substep(float const dt, int const frame);

    // Make the velocity on the grid to be divergence free
    void pressureProjection(bool print);

    void updateEmitter();

    void applyDirichletVelocity(Vec3SGrid& vecGrid, int frame);

    void swap();

    void addGravity(float const dt);

    void advectDensity(float const dt);

    void advectVelocity(float const dt, int const frame);

    void writeVDBs(int const frame);
    struct BoundaryOp {
        BoundaryOp(Vec3SGrid::ConstPtr dirichletVelocity,
                   FloatGrid::ConstPtr dirichletPressure,
                   float const voxelSize) :
                   dirichletVelocity(dirichletVelocity),
                   dirichletPressure(dirichletPressure),
                   voxelSize(voxelSize) {}

        void operator()(const openvdb::Coord& ijk,
                        const openvdb::Coord& neighbor,
                        double& source,
                        double& diagonal) const
        {
            float const dirichletBC = 0.f;
            bool isNeumannPressure = dirichletVelocity->tree().isValueOn(neighbor);
            auto vNgbr = dirichletVelocity->tree().getValue(neighbor);
            // bool isDirichletPressure = dirichletPressure->tree().isValueOn(neighbor);

            // TODO: Double check this:
            if (isNeumannPressure) {
                double delta = 0.0;
                // Neumann pressure from bbox
                if (neighbor.x() + 1 == ijk.x() /* left x-face */) {
                    delta += vNgbr[0];
                }
                if (neighbor.x() - 1 == ijk.x() /* right x-face */) {
                    delta -= vNgbr[0];
                }
                if (neighbor.y() + 1 == ijk.y() /* bottom y-face */) {
                    delta += vNgbr[1];
                }
                if (neighbor.y() - 1 == ijk.y() /* top y-face */) {
                    delta -= vNgbr[1];
                }
                if (neighbor.z() + 1 == ijk.z() /* back z-face */) {
                    delta +=  vNgbr[2];
                }
                if (neighbor.z() - 1 == ijk.z() /* front z-face */) {
                    delta -=  vNgbr[2];
                }
                // Note: in the SOP_OpenVDB_Remove_Divergence, we need to multiply
                // this by 0.5, because the gradient that's used is using
                // central-differences in a collocated grid, instead of the staggered one.
                source += delta / voxelSize;
            } else {
                // Dirichlet pressure
                if (neighbor.x() + 1 == ijk.x() /* left x-face */) {
                    diagonal -= 1.0;
                    source -= dirichletBC;
                }
                if (neighbor.x() - 1 == ijk.x() /* right x-face */) {
                    diagonal -= 1.0;
                    source -= dirichletBC;
                }
                if (neighbor.y() + 1 == ijk.y() /* bottom y-face */) {
                    diagonal -= 1.0;
                    source -= dirichletBC;
                }
                if (neighbor.y() - 1 == ijk.y() /* top y-face */) {
                    diagonal -= 1.0;
                    source -= dirichletBC;
                }
                if (neighbor.z() + 1 == ijk.z() /* back z-face */) {
                    diagonal -= 1.0;
                    source -= dirichletBC;
                }
                if (neighbor.z() - 1 == ijk.z() /* front z-face */) {
                    diagonal -= 1.0;
                    source -= dirichletBC;
                }
            }
        }

        Vec3SGrid::ConstPtr dirichletVelocity;
        FloatGrid::ConstPtr dirichletPressure;
        float voxelSize;
    };


    // Apply Gravity Functor. Meant to be used with
    // foreach in LeafManager
    struct ApplyGravityOp
    {
        ApplyGravityOp(float const dt, Vec3s const gravity) : dt(dt), gravity(gravity) {}

        template <typename T>
        void operator()(T &leaf, size_t) const
        {
            for (typename T::ValueOnIter iter = leaf.beginValueOn(); iter; ++iter) {
                Vec3s newVal = *iter  + dt * gravity;
                iter.setValue(newVal);
            }
        }

        Vec3s const gravity;
        float const dt;
    };// ApplyGravityOp


    struct ApplyDirichletVelocityOp
    {
        ApplyDirichletVelocityOp(Vec3SGrid::Ptr dirichletVelocity) :
            dirichletVelocity(dirichletVelocity) {}

        template <typename T>
        void operator()(T &leaf, size_t) const
        {
            auto drcVelAcc = dirichletVelocity->getConstAccessor();
            for (typename T::ValueOnIter iter = leaf.beginValueOn(); iter; ++iter) {
                auto ijk = iter.getCoord();
                if (drcVelAcc.isValueOn(ijk)) {
                    iter.setValue(drcVelAcc.getValue(ijk));
                }
            }
        }

        Vec3SGrid::ConstPtr dirichletVelocity;
    };// ApplyDirichletVelocityOp


    // Update the emitter by taking the maximum density with the current density.
    struct UpdateEmitterOp
    {
        UpdateEmitterOp(FloatGrid::Ptr emitter, FloatGrid::Ptr density) : emitter(emitter), density(density) {}

        template <typename T>
        void operator()(T &leaf, size_t) const
        {
            auto densityAcc = density->getAccessor();
            for (typename T::ValueOnIter iter = leaf.beginValueOn(); iter; ++iter) {
                auto ijk = iter.getCoord();
                if (densityAcc.isValueOn(ijk)) {
                    float newVal = std::max(*iter, densityAcc.getValue(ijk));
                    densityAcc.setValue(ijk, newVal);
                }
            }
        }

        FloatGrid::Ptr emitter;
        FloatGrid::Ptr density; // current density
    };// UpdateEmitterOp


    float mVoxelSize = 0.1f;
    Vec3s mGravity = Vec3s(0.f, -9.8f, 0.f);
    math::Transform::Ptr mXform;

    FloatGrid::Ptr mDirichletPressure;

    FloatGrid::Ptr mDivBefore;
    FloatGrid::Ptr mDivAfter;

    FloatGrid::Ptr mEmitter;
    Vec3SGrid::Ptr mDirichletVelocity;

    FloatGrid::Ptr mDensityCurr;
    FloatGrid::Ptr mDensityNext;
    Vec3SGrid::Ptr mVCurr;
    Vec3SGrid::Ptr mVNext;
    FloatGrid::Ptr mPressure;
    BoolGrid::Ptr mDomainMaskGrid;
};


SmokeSolver::SmokeSolver(float const voxelSize) : mVoxelSize(voxelSize)
{
    initialize();
}


void
SmokeSolver::swap()
{
    std::cout << "Swap" << std::endl;
    mDensityCurr = mDensityNext;
    mDensityCurr->setName("density_curr");
    mVCurr = mVNext;
    mVCurr->setName("vel_curr");
}


void
SmokeSolver::applyDirichletVelocity(Vec3SGrid& vecGrid, int frame)
{
    std::cout << "apply dirichlet velocity begins" << std::endl;
    auto vecAcc = vecGrid.getAccessor();
    // mVNext = mVCurr->deepCopy();
    //auto vNextAcc = mVNext->getAccessor();
    auto drcAcc = mDirichletVelocity->getConstAccessor();
    bool print = true;
    for (auto iter = vecGrid.beginValueOn(); iter; ++iter) {
        auto ijk = iter.getCoord();
        if (drcAcc.isValueOn(ijk)) {
           // vecAcc.setValue(ijk, *iter);
           //if (print) {
           //std::cout << "\t\tdebug apply dirichlet velocity frame = " << frame << "\tijk = " << ijk << "\t *drc val = " << drcAcc.getValue(ijk) << "\tbefore val = " << vNextAcc.getValue(ijk) << "\n";
           //}
           print = false;
           //vNextAcc.setValue(ijk, drcAcc.getValue(ijk));
           vecAcc.setValue(ijk, drcAcc.getValue(ijk));
        }
    }

    //std::ostringstream ostr;
    //ostr << "debug_apply_dirichet_velocity" << "_" << frame << ".vdb";
    //std::cerr << "\tWriting " << ostr.str() << std::endl;
    //openvdb::io::File file(ostr.str());
    //openvdb::GridPtrVec grids;
    //grids.push_back(mVCurr);
    //grids.push_back(mVNext);
    //mVCurr->setName("vel_curr");
    //mVCurr->setName("vel_next");
    //file.write(grids);
    //file.close();
    // tree::LeafManager<Vec3STree> lm(vecGrid.tree());
    // SmokeSolver::ApplyDirichletVelocityOp op(mDirichletVelocity);
    // lm.foreach(op);
}

void
SmokeSolver::updateEmitter()
{
    std::cout << "update emitter" << std::endl;
    tree::LeafManager<FloatTree> lm(mEmitter->tree());
    UpdateEmitterOp op(mEmitter, mDensityCurr);
    lm.foreach(op);
    std::cout << "done with update emitter" << std::endl;
}


void
SmokeSolver::initialize() {
    using BBox = math::BBox<Vec3s>;

    mXform = math::Transform::createLinearTransform(mVoxelSize);
    float const padding = 4.f * mVoxelSize;
    float const centerY = 3.f;

    // Create an emitter and an emitter velocity
    auto minEmtW = Vec3s(0.f, 2.5f, 2.5f);
    auto maxEmtW = Vec3s(2 * mVoxelSize, 3.5f, 3.5f);
    Coord minEmtCoord = mXform->worldToIndexNodeCentered(minEmtW);
    Coord maxEmtCoord = mXform->worldToIndexNodeCentered(maxEmtW);
    mEmitter = FloatGrid::create(/*bg = */0.f);
    mEmitter->denseFill(CoordBBox(minEmtCoord, maxEmtCoord), /* value = */ 2.0, /*active = */ true);
    mEmitter->setTransform(mXform);
    mEmitter->setName("emitter");

    auto minDVLft = Vec3s(-padding, -padding, -padding);
    auto maxDVLft = Vec3s(2 * mVoxelSize, 6.f, 6.f);
    Coord minDVLftCoord = mXform->worldToIndexNodeCentered(minDVLft);
    Coord maxDVLftCoord = mXform->worldToIndexNodeCentered(maxDVLft);
    mDirichletVelocity = Vec3SGrid::create(/* bg = */ Vec3s(0.f, 0.f, 0.f));
    mDirichletVelocity->denseFill(CoordBBox(minDVLftCoord, maxDVLftCoord), /* value = */ Vec3s(1.f, 0.f, 0.f), /*active = */ true);
    mDirichletVelocity->setGridClass(GRID_STAGGERED);
    mDirichletVelocity->setTransform(mXform);
    mDirichletVelocity->setName("dirichlet_velocity");

    // Create interior mask
    auto minBBox = Vec3s(0.f, 0.f, 0.f);
    auto maxBBox = Vec3s(14.f, 6.f, 6.f);
    Coord minBBoxIntrCoord = mXform->worldToIndexNodeCentered(minBBox);
    Coord maxBBoxIntrCoord = mXform->worldToIndexNodeCentered(maxBBox);
    mDomainMaskGrid = BoolGrid::create(false /* background */);
    mDomainMaskGrid->denseFill(CoordBBox(minBBoxIntrCoord, maxBBoxIntrCoord), /* value = */ true, /* active = */ true);
    mDomainMaskGrid->setTransform(mXform);
    mDomainMaskGrid->setName("domain_mask");

    // Set up density and velocity grid. Need to take the collider out.
    mDensityCurr = FloatGrid::create(/*bg = */0.f);
    mDensityCurr->denseFill(CoordBBox(minBBoxIntrCoord, maxBBoxIntrCoord), /* value = */ 0.f, /* active = */ true);
    mDensityCurr->setTransform(mXform);
    mDensityCurr->setName("density_curr");
    mDensityCurr->topologyUnion(*mDomainMaskGrid);
    mDensityCurr->topologyIntersection(*mDomainMaskGrid);

    mVCurr = Vec3SGrid::create(/*bg = */Vec3s(0.f, 0.f, 0.f));
    mVCurr->denseFill(CoordBBox(minBBoxIntrCoord, maxBBoxIntrCoord), /* value = */ Vec3s(1.f, 0.f, 0.f), /* active = */ true);
    mVCurr->setTransform(mXform);
    mVCurr->setName("vel_curr");
    mVCurr->setGridClass(GRID_STAGGERED);
    mVCurr->topologyUnion(*mDomainMaskGrid);
    mVCurr->topologyIntersection(*mDomainMaskGrid);
    mVCurr->tree().voxelizeActiveTiles();

    updateEmitter();
    applyDirichletVelocity(*mVCurr, -1);
}


void
SmokeSolver::addGravity(float const dt) {
    tree::LeafManager<Vec3STree> lm(mVCurr->tree());
    SmokeSolver::ApplyGravityOp op(dt, mGravity);
    lm.foreach(op);
}


void
SmokeSolver::advectDensity(float const dt)
{
    std::cout << "Advect Density" << std::endl;
    using MaskGridType = BoolGrid;
    using AdvT = openvdb::tools::VolumeAdvection<Vec3fGrid, true /* Staggered */>;
    using SamplerT = openvdb::tools::Sampler<1>;

    AdvT advection(*mVCurr);
    advection.setIntegrator(tools::Scheme::MAC);
    advection.setLimiter(tools::Scheme::REVERT);
    advection.setSubSteps(1);

    mDensityNext = advection.advect<FloatGrid, BoolGrid, SamplerT>(*mDensityCurr, *mDomainMaskGrid, dt);
    mDensityNext->setName("density_next");

    // Debug
    // {//test advect with a mask
    //     Vec3fGrid velocity(Vec3f(1.0f, 0.0f, 0.0f));
    //     FloatGrid::Ptr density0 = FloatGrid::create(0.0f), density1;
    //     density0->fill(CoordBBox(Coord(0),Coord(6)), 1.0f);

    //     BoolGrid::Ptr mask = BoolGrid::create(false);
    //     mask->fill(CoordBBox(Coord(4,0,0),Coord(30,8,8)), true);

    //     AdvT a(*mVCurr);
    //     a.setGrainSize(0);
    //     a.setIntegrator(tools::Scheme::MAC);
    //     //a.setIntegrator(tools::Scheme::BFECC);
    //     //a.setIntegrator(tools::Scheme::RK4);
    //     for (int i=1; i<=240; ++i) {
    //         // density1 = a.advect<FloatGrid, BoolGrid, SamplerT>(*density0, *mask, 0.1f);
    //         mDensityNext = a.advect<FloatGrid, BoolGrid, SamplerT>(*mDensityCurr, *mDomainMaskGrid, 0.1f);
    //         std::ostringstream ostr;
    //         ostr << "debugAdvectMask" << "_" << i << ".vdb";
    //         std::cerr << "Writing " << ostr.str() << std::endl;
    //         openvdb::io::File file(ostr.str());
    //         openvdb::GridPtrVec grids;
    //         // grids.push_back(density1);
    //         grids.push_back(mDensityNext);
    //         file.write(grids);
    //         // density0 = density1;
    //         mDensityCurr = mDensityNext;
    //     }
    // }
    // std::ostringstream ostr;
    // ostr << "densityAdvectMask" << ".vdb";
    // std::cerr << "Writing " << ostr.str() << std::endl;
    // openvdb::io::File file(ostr.str());
    // openvdb::GridPtrVec grids;
    // grids.push_back(mDensityNext);
    // file.write(grids);
    // file.close();
}

void
SmokeSolver::advectVelocity(float const dt, const int frame)
{
    std::cout << "Advect Velocity" << std::endl;
    using AdvT = openvdb::tools::VolumeAdvection<Vec3SGrid, true /* staggered */>;
    using SamplerT = openvdb::tools::Sampler<1, true /* staggered */>;

    AdvT advection(*mVCurr);
    advection.setIntegrator(tools::Scheme::MAC);
    advection.setLimiter(tools::Scheme::REVERT);
    advection.setSubSteps(1);

    mVNext = advection.advect<Vec3SGrid, BoolGrid, SamplerT>(*mVCurr, *mDomainMaskGrid, dt);
    // mVNext = advection.advect<Vec3SGrid, SamplerT>(*mVCurr, dt);
    mVNext->setGridClass(GRID_STAGGERED);
    mVNext->setName("vel_next");

    std::ostringstream ostr;
    ostr << "debug_velocity_advection" << "_" << frame << ".vdb";
    std::cerr << "\tvelocity advection::Writing " << ostr.str() << std::endl;
    std::cout << "\tvelocity advection::mVCurr->background() = " << mVCurr->background() << std::endl;
    openvdb::io::File file(ostr.str());
    openvdb::GridPtrVec grids;
    grids.push_back(mVNext);
    grids.push_back(mVCurr);
    grids.push_back(mDomainMaskGrid);
    file.write(grids);
    file.close();
}


void
SmokeSolver::pressureProjection(bool print) {
    using TreeType = FloatTree;
    using ValueType = TreeType::ValueType;
    using MaskGridType = BoolGrid;
    using PCT = openvdb::math::pcg::JacobiPreconditioner<openvdb::tools::poisson::LaplacianMatrix>;

    ValueType const zero = zeroVal<ValueType>();
    double const epsilon = math::Delta<ValueType>::value();

    mDivBefore = tools::divergence(*mVCurr);
    mDivBefore->setName("div_before");

    MaskGridType* domainMaskGrid = new MaskGridType(*mDivBefore); // match input grid's topology
    tools::erodeActiveValues(domainMaskGrid->tree(), /*iterations=*/1, tools::NN_FACE, tools::IGNORE_TILES);
    // domainMaskGrid->topologyDifference(*mDirichletPressure);
    domainMaskGrid->topologyDifference(*mDirichletVelocity);

    math::pcg::State state = math::pcg::terminationDefaults<ValueType>();
    state.iterations = 100000;
    state.relativeError = state.absoluteError = epsilon;
    SmokeSolver::BoundaryOp bop(mDirichletVelocity, mDirichletPressure, mVoxelSize);
    util::NullInterrupter interrupter;
    FloatTree::Ptr fluidPressure = tools::poisson::solveWithBoundaryConditionsAndPreconditioner<PCT>(
        mDivBefore->tree(), domainMaskGrid->tree(), bop, state, interrupter, /*staggered=*/true);
    FloatGrid::Ptr fluidPressureGrid = FloatGrid::create(fluidPressure);
    // Note: need to dilate in order to do one-sided difference
    // because we use a staggered grid velocity field.
    tools::dilateActiveValues(*fluidPressure, /*iterations=*/1, tools::NN_FACE, tools::IGNORE_TILES);

    fluidPressureGrid->setTransform(mXform);
    mPressure = fluidPressureGrid->copy();
    mPressure->setName("pressure");

    auto vCurrAcc = mVCurr->getAccessor();
    // auto vNextAcc = mVNext->getAccessor();
    auto pressureAcc = fluidPressureGrid->getAccessor();
    // Note: I'm modifying vCurr
    for (auto iter = mVCurr->beginValueOn(); iter; ++iter) {
        math::Coord ijk = iter.getCoord();
        Vec3s gradijk;
        gradijk[0] = pressureAcc.getValue(ijk) - pressureAcc.getValue(ijk.offsetBy(-1, 0, 0));
        gradijk[1] = pressureAcc.getValue(ijk) - pressureAcc.getValue(ijk.offsetBy(0, -1, 0));
        gradijk[2] = pressureAcc.getValue(ijk) - pressureAcc.getValue(ijk.offsetBy(0, 0, -1));

        // This is only multiplied by mVoxelSize because in the computation of gradijk, I don't divide by mVoxelSize.
        auto val = vCurrAcc.getValue(ijk) - gradijk * mVoxelSize;
        vCurrAcc.setValue(ijk, val);
    }

    std::cout << "Projection Success: " << state.success << "\n";
    std::cout << "Iterations: " << state.iterations << "\n";
    std::cout << "Relative error: " << state.relativeError << "\n";
    std::cout << "Absolute error: " << state.absoluteError << "\n";
}


void
SmokeSolver::substep(float const dt, int const frame) {
    mVCurr->setName("vel_curr");
    // updateEmitter();
    // addGravity(dt);
    applyDirichletVelocity(*mVCurr, frame);
    //pressureProjection(true /* print */);
    //applyDirichletVelocity(*mVCurr);
    //advectDensity(dt);
    advectVelocity(dt, frame);
    applyDirichletVelocity(*mVNext, frame);
}


void
SmokeSolver::foobar() {
    float const dt = 1.f/24.f;
    for (int frame = 0; frame < 10; ++frame) {
        std::cout << "\n====== foobar frame " << frame << " ======" << std::endl;
        // updateEmitter();
        {
            std::ostringstream ostr;
            ostr << "before advect density" << "_" << frame << ".vdb";
            openvdb::io::File file(ostr.str());
            openvdb::GridPtrVec grids;
            grids.push_back(mVCurr);
            file.write(grids);
            file.close();
        }
        advectDensity(dt);
        {
            std::ostringstream ostr;
            ostr << "before advect velocity" << "_" << frame << ".vdb";
            openvdb::io::File file(ostr.str());
            openvdb::GridPtrVec grids;
            grids.push_back(mVCurr);
            file.write(grids);
            file.close();
        }
        advectVelocity(dt, frame);
        swap();
        applyDirichletVelocity(*mVCurr, frame);
        writeVDBs(frame);
    }
}


void
SmokeSolver::render() {
    float const dt = 1.f/24.f;
    for (int frame = 0; frame < 10; ++frame) {
        std::cout << "\nframe = " << frame << "\n";
        int const numSubstep = 4;
        for (int i = 0; i < numSubstep; ++i) {
            std::cout << "\tsubstep = " << i << std::endl;
            substep(dt / numSubstep, frame);
            //if (i < numSubstep - 1) {
            swap();
            //}
        }
        writeVDBs(frame);
    }
}




void
SmokeSolver::writeVDBs(int const frame) {
    std::cout << "Write VDBs" << std::endl;
    std::ostringstream ss;
    ss << "smoke_" << std::setw(3) << std::setfill('0') << frame << ".vdb";
    std::string fileName(ss.str());
    io::File file(fileName.c_str());

    openvdb::GridPtrVec grids;
    grids.push_back(mEmitter);
    grids.push_back(mDirichletVelocity);
    grids.push_back(mDensityCurr);
    grids.push_back(mDensityNext);
    grids.push_back(mVCurr);
    // grids.push_back(mVNext);
    // grids.push_back(mDirichletPressure);
    grids.push_back(mDivBefore);
    grids.push_back(mDivAfter);
    grids.push_back(mPressure);
    grids.push_back(mDomainMaskGrid);

    file.write(grids);
    file.close();
}
}