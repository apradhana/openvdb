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
#include <openvdb/points/PointAdvect.h> // for advectPoints
#include <openvdb/points/PointAttribute.h> // for appendAttribute
#include <openvdb/points/PointDataGrid.h> // for PointDataGrid
#include <openvdb/points/PointRasterizeTrilinear.h> // for rasterizing to the grid
#include <openvdb/points/PointSample.h> // for PointSample
#include <openvdb/points/PointScatter.h> // for point sampling
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

private:

    void initialize();

    void substep(float const dt);

    // Make the velocity on the grid to be divergence free
    void pressureProjection(bool print);

    void updateEmitter();

    void gridVelocityUpdate(float const dt);

    void velocityBCCorrection(Vec3SGrid& vecGrid);

    void addGravity(float const dt);

    void writeVDBs(int const frame);
    struct BoundaryOp {
        BoundaryOp(float const voxelSize,
                   FloatGrid::Ptr collider,
                   Vec3SGrid::Ptr vCurr) :
            voxelSize(voxelSize),
            collider(collider),
            vCurr(vCurr) {}

        void operator()(const openvdb::Coord& ijk,
                        const openvdb::Coord& neighbor,
                        double& source,
                        double& diagonal) const
        {
            float const dirichletBC = 0.f;
            bool isInsideCollider = collider->tree().isValueOn(neighbor);
            auto vNgbr = vCurr->tree().getValue(neighbor);

            // TODO: Double check this:
            if (isInsideCollider) {
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

        float voxelSize;
        FloatGrid::Ptr collider;
        Vec3SGrid::Ptr vCurr;
    };


    struct FlipUpdateOp
    {
        explicit FlipUpdateOp(Index64 const velAtrIdx,
                              Index64 const vPicAtrIdx,
                              Index64 const vFlipAtrIdx,
                              float const alpha)
                              : velAtrIdx(velAtrIdx),
                                vPicAtrIdx(vPicAtrIdx),
                                vFlipAtrIdx(vFlipAtrIdx),
                                alpha(alpha) { }

        void operator()(const tree::LeafManager<points::PointDataTree>::LeafRange& range) const {
            for (auto leafIter = range.begin(); leafIter; ++leafIter) {
                points::AttributeArray& velArray = leafIter->attributeArray(velAtrIdx);
                points::AttributeArray const& vPicArray = leafIter->constAttributeArray(vPicAtrIdx);
                points::AttributeArray const& vFlipArray = leafIter->constAttributeArray(vFlipAtrIdx);
                points::AttributeWriteHandle<Vec3s> velHandle(velArray);
                points::AttributeHandle<Vec3s> vPicHandle(vPicArray);
                points::AttributeHandle<Vec3s> vFlipHandle(vFlipArray);
                // Iterate over active indices in the leaf.
                for (auto indexIter = leafIter->beginIndexOn(); indexIter; ++indexIter) {
                    auto vPic = vPicHandle.get(*indexIter);
                    auto vFlip = vFlipHandle.get(*indexIter);
                    auto newVel = alpha * (vPic + vFlip) + (1 - alpha) * vPic;
                    velHandle.set(*indexIter, newVel);
                }
            }
        }

        Index64 velAtrIdx;
        Index64 vPicAtrIdx;
        Index64 vFlipAtrIdx;
        float alpha;
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
                float newVal = std::max(*iter, densityAcc.getValue(ijk));
                densityAcc.setValue(ijk, newVal);
            }
        }

        FloatGrid::Ptr emitter;
        FloatGrid::Ptr density; // current density
    };// UpdateEmitterOp


    // Compute the difference between vNext and the original rasterized
    // vCurr (before the addition of gravity). To be used with foreach in LeafManager.
    struct ComputeFlipVelocityOp
    {
        ComputeFlipVelocityOp(Vec3SGrid::Ptr vCurr,
                              Vec3SGrid::Ptr vNext,
                              float const dt,
                              Vec3s const gravity) :
                              vCurr(vCurr),
                              vNext(vNext),
                              dt(dt),
                              gravity(gravity) {}

        template <typename T>
        void operator()(T &leaf, size_t) const
        {
            auto vCurrAcc = vCurr->getAccessor();
            auto vNextAcc = vNext->getAccessor();
            for (typename T::ValueOnIter iter = leaf.beginValueOn(); iter; ++iter) {
                auto ijk = iter.getCoord();
                Vec3s val = vNextAcc.getValue(ijk) - vCurrAcc.getValue(ijk) - dt * gravity;
                iter.setValue(val);
            }
        }

        Vec3SGrid::Ptr vCurr;
        Vec3SGrid::Ptr vNext;
        Vec3s const gravity;
        float const dt;
    };// ComputeFlipVelocityOp


    float mVoxelSize = 0.1f;
    Vec3s mGravity = Vec3s(0.f, -9.8f, 0.f);
    math::Transform::Ptr mXform;

    FloatGrid::Ptr mCollider;

    FloatGrid::Ptr mDivBefore;
    FloatGrid::Ptr mDivAfter;

    FloatGrid::Ptr mEmitter;
    Vec3SGrid::Ptr mDirichletVelocity;

    FloatGrid::Ptr mDensityCurr;
    FloatGrid::Ptr mDensityNext;
    Vec3SGrid::Ptr mVCurr;
    Vec3SGrid::Ptr mVNext;
    FloatGrid::Ptr mPressure;
};


SmokeSolver::SmokeSolver(float const voxelSize) : mVoxelSize(voxelSize)
{
    initialize();
}

void
SmokeSolver::updateEmitter()
{
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
    auto maxEmtW = Vec3s(3.f, 3.5f, 3.5f);
    Coord minEmtCoord = mXform->worldToIndexNodeCentered(minEmtW);
    Coord maxEmtCoord = mXform->worldToIndexNodeCentered(maxEmtW);
    mEmitter = FloatGrid::create(/*bg = */0.f);
    mEmitter->denseFill(CoordBBox(minEmtCoord, maxEmtCoord), /* value = */ 2.0, /*active = */ true);
    mEmitter->setTransform(mXform);
    mEmitter->setName("emitter");

    // Create Dirichlet Velocity (Neumann-pressure)
    auto minDirichletVelW = Vec3s(0.f, 0.f, 0.f);
    auto maxDirichletVelW = Vec3s(2 * mVoxelSize, 6.f, 6.f);
    Coord minDirichletVelCoord = mXform->worldToIndexNodeCentered(minDirichletVelW);
    Coord maxDirichletVelCoord = mXform->worldToIndexNodeCentered(maxDirichletVelW);
    mDirichletVelocity = Vec3SGrid::create(/* bg = */ Vec3s(0.f, 0.f, 0.f));
    mDirichletVelocity->denseFill(CoordBBox(minDirichletVelCoord, maxDirichletVelCoord), /* value = */ Vec3s(1.f, 0.f, 0.f), /*active = */ true);
    mDirichletVelocity->setTransform(mXform);
    mDirichletVelocity->setName("dirichlet_velocity");


    // Create collider. Collider is bounding box minus interior plus sphere.
    auto minBBox = Vec3s(0.f, 0.f, 0.f);
    auto maxBBox = Vec3s(14.f + mVoxelSize, 6.f + mVoxelSize, 6.f + mVoxelSize);
    Coord minBBoxIntrCoord = mXform->worldToIndexNodeCentered(minBBox);
    Coord maxBBoxIntrCoord = mXform->worldToIndexNodeCentered(maxBBox);

    FloatGrid::Ptr negativeSpace = FloatGrid::create(/*bg = */0.f);
    negativeSpace->denseFill(CoordBBox(minBBoxIntrCoord, maxBBoxIntrCoord), /*value = */ 1.0, /*active = */ true);
    negativeSpace->setTransform(mXform);

    Vec3s const center(5.f, 3.f, 3.f);
    float const radius = 1.f;
    FloatGrid::Ptr sphere = tools::createLevelSetSphere<FloatGrid>(radius, center, mVoxelSize);
    sphere->setTransform(mXform);
    sphere->setName("sphere");

    Vec3s minBBoxvec = Vec3s(-padding, -padding, -padding);
    Vec3s maxBBoxvec = Vec3s(14.f + padding + mVoxelSize, 6.f + padding + mVoxelSize, 6.f + padding + mVoxelSize);
    Coord minBBoxcoord = mXform->worldToIndexNodeCentered(minBBoxvec);
    Coord maxBBoxcoord = mXform->worldToIndexNodeCentered(maxBBoxvec);
    mCollider = FloatGrid::create(/*bg = */0.f);
    mCollider->denseFill(CoordBBox(minBBoxcoord, maxBBoxcoord), /*value = */ 1.0, /*active = */ true);
    mCollider->topologyDifference(*negativeSpace);
    mCollider->topologyUnion(*sphere);
    mCollider->setTransform(mXform);
    mCollider->setName("collider");


    // Set up density and velocity grid. Need to take the collider out.
    mDensityCurr = FloatGrid::create(/*bg = */0.f);
    mDensityCurr->denseFill(CoordBBox(minBBoxIntrCoord, maxBBoxIntrCoord), /* value = */ 0.f, /* active = */ true);
    mDensityCurr->setTransform(mXform);
    mDensityCurr->setName("density_curr");

    mDensityNext = FloatGrid::create(/*bg = */0.f);
    mDensityNext->denseFill(CoordBBox(minBBoxIntrCoord, maxBBoxIntrCoord), /* value = */ 0.f, /* active = */ true);
    mDensityNext->setTransform(mXform);
    mDensityNext->setName("density_next");

    mVCurr = Vec3SGrid::create(/*bg = */Vec3s(0.f, 0.f, 0.f));
    mVCurr->denseFill(CoordBBox(minBBoxIntrCoord, maxBBoxIntrCoord), /* value = */ Vec3s(0.f, 0.f, 0.f), /* active = */ true);
    mVCurr->setTransform(mXform);
    mVCurr->setName("vel_curr");

    mVNext = Vec3SGrid::create(/*bg = */Vec3s(0.f, 0.f, 0.f));
    mVNext->denseFill(CoordBBox(minBBoxIntrCoord, maxBBoxIntrCoord), /* value = */ Vec3s(0.f, 0.f, 0.f), /* active = */ true);
    mVNext->setTransform(mXform);
    mVNext->setName("vel_next");
}


void
SmokeSolver::addGravity(float const dt) {
    tree::LeafManager<Vec3STree> lm(mVCurr->tree());
    SmokeSolver::ApplyGravityOp op(dt, mGravity);
    lm.foreach(op);
}


void
SmokeSolver::velocityBCCorrection(Vec3SGrid& vecGrid) {
    auto acc = vecGrid.getAccessor();
    auto cldrAcc = mCollider->getAccessor();

    for (auto iter = vecGrid.beginValueOn(); iter; ++iter) {
        math::Coord ijk = iter.getCoord();
        math::Coord im1jk = ijk.offsetBy(-1, 0, 0);
        math::Coord ip1jk = ijk.offsetBy(1, 0, 0);
        math::Coord ijm1k = ijk.offsetBy(0, -1, 0);
        math::Coord ijp1k = ijk.offsetBy(0, 1, 0);
        math::Coord ijkm1 = ijk.offsetBy(0, 0, -1);
        math::Coord ijkp1 = ijk.offsetBy(0, 0, 1);

        if (cldrAcc.isValueOn(im1jk) || cldrAcc.isValueOn(ip1jk)) {
            auto val = acc.getValue(ijk);
            Vec3s newVal = Vec3s(0, val[1], val[2]);
            acc.setValue(ijk, newVal);
        }
        if (cldrAcc.isValueOn(ijm1k) || cldrAcc.isValueOn(ijp1k)) {
            auto val = acc.getValue(ijk);
            Vec3s newVal = Vec3s(val[0], 0, val[2]);
            acc.setValue(ijk, newVal);
        }
        if (cldrAcc.isValueOn(ijkm1) || cldrAcc.isValueOn(ijkp1)) {
            auto val = acc.getValue(ijk);
            Vec3s newVal = Vec3s(val[0], val[1], 0);
            acc.setValue(ijk, newVal);
        }
    }
}


void
SmokeSolver::pressureProjection(bool print) {
    using TreeType = FloatTree;
    using ValueType = TreeType::ValueType;
    using MaskGridType = BoolGrid;
    using PCT = openvdb::math::pcg::JacobiPreconditioner<openvdb::tools::poisson::LaplacianMatrix>;

    ValueType const zero = zeroVal<ValueType>();
    double const epsilon = math::Delta<ValueType>::value();

    BoolTree::Ptr interiorMask(new BoolTree(false));
    interiorMask->topologyUnion(mVCurr->tree());
    tools::erodeActiveValues(*interiorMask, /*iterations=*/1, tools::NN_FACE, tools::IGNORE_TILES);
    BoolGrid::Ptr interiorGrid = BoolGrid::create(interiorMask);
    interiorGrid->setTransform(mXform);
    // mInterior = interiorGrid->copy();
    // mInterior->setName("interior");

    mDivBefore = tools::divergence(*mVCurr);
    mDivBefore->setName("div_before");

    MaskGridType* domainMaskGrid = new MaskGridType(*mDivBefore); // match input grid's topology
    tools::erodeActiveValues(domainMaskGrid->tree(), /*iterations=*/1, tools::NN_FACE, tools::IGNORE_TILES);
    domainMaskGrid->topologyDifference(*mCollider);

    math::pcg::State state = math::pcg::terminationDefaults<ValueType>();
    state.iterations = 100000;
    state.relativeError = state.absoluteError = epsilon;
    SmokeSolver::BoundaryOp bop(mVoxelSize, mCollider, mVCurr);
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
    auto vNextAcc = mVNext->getAccessor();
    auto pressureAcc = fluidPressureGrid->getAccessor();
    for (auto iter = mVCurr->beginValueOn(); iter; ++iter) {
        math::Coord ijk = iter.getCoord();
        Vec3s gradijk;
        gradijk[0] = pressureAcc.getValue(ijk) - pressureAcc.getValue(ijk.offsetBy(-1, 0, 0));
        gradijk[1] = pressureAcc.getValue(ijk) - pressureAcc.getValue(ijk.offsetBy(0, -1, 0));
        gradijk[2] = pressureAcc.getValue(ijk) - pressureAcc.getValue(ijk.offsetBy(0, 0, -1));

        // This is only multiplied by mVoxelSize because in the computation of gradijk, I don't divide by mVoxelSize.
        auto val = vCurrAcc.getValue(ijk) - gradijk * mVoxelSize;
        vNextAcc.setValue(ijk, val);
    }

    std::cout << "Projection Success: " << state.success << "\n";
    std::cout << "Iterations: " << state.iterations << "\n";
    std::cout << "Relative error: " << state.relativeError << "\n";
    std::cout << "Absolute error: " << state.absoluteError << "\n";
}


void
SmokeSolver::gridVelocityUpdate(float const dt) {
    addGravity(dt);
    velocityBCCorrection(*mVCurr);
    pressureProjection(false /* print */);
    velocityBCCorrection(*mVNext);
}


void
SmokeSolver::substep(float const dt) {
    updateEmitter();
    // gridVelocityUpdate(dt);
}


void
SmokeSolver::render() {
    float const dt = 1.f/24.f;
    for (int frame = 0; frame < 100; ++frame) {
        std::cout << "\nframe = " << frame << "\n";
        substep(dt);
        writeVDBs(frame);
    }
}


void
SmokeSolver::writeVDBs(int const frame) {
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
    grids.push_back(mVNext);
    grids.push_back(mCollider);
    grids.push_back(mDivBefore);
    grids.push_back(mDivAfter);
    grids.push_back(mPressure);

    file.write(grids);
    file.close();
}
}