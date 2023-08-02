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

    // Rasterize particle velocity to the grid
    void particlesToGrid();

    // FLIP update: Interpolate the delta of velocity update (v_np1 - v_n)
    // back to the particle
    void gridToParticles();
    void updateParticles(float const dt);
    void updateParticlesVelocity();

    // Update particle position based on velocity on the grid
    void advectParticles(float const dt);

    // Make the velocity on the grid to be divergence free
    void pressureProjection(bool print);

    void gridVelocityUpdate(float const dt);

    void velocityBCCorrection(Vec3SGrid& vecGrid);

    void addGravity(float const dt);
    void computeFlipVelocity(float const dt);

    void writeVDBs(int const frame);
    void writeVDBsVerbose(int const frame);

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
                    delta += /* voxelSize * */ vNgbr[0];
                }
                if (neighbor.x() - 1 == ijk.x() /* right x-face */) {
                    delta -= /* voxelSize * */ vNgbr[0];
                }
                if (neighbor.y() + 1 == ijk.y() /* bottom y-face */) {
                    delta += /* voxelSize * */ vNgbr[1];
                }
                if (neighbor.y() - 1 == ijk.y() /* top y-face */) {
                    delta -= /* voxelSize * */ vNgbr[1];
                }
                if (neighbor.z() + 1 == ijk.z() /* back z-face */) {
                    delta += /* voxelSize *  */ vNgbr[2];
                }
                if (neighbor.z() - 1 == ijk.z() /* front z-face */) {
                    delta -= /* voxelSize *  */ vNgbr[2];
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
    int mPointsPerVoxel = 8;
    math::Transform::Ptr mXform;

    points::PointDataGrid::Ptr mPoints;
    FloatGrid::Ptr mBBoxLS;
    FloatGrid::Ptr mCollider;

    FloatGrid::Ptr mDivBefore;
    FloatGrid::Ptr mDivAfter;

    FloatGrid::Ptr mDensity;
    Vec3SGrid::Ptr mVCurr;
    Vec3SGrid::Ptr mVNext;
    Vec3SGrid::Ptr mVDiff; // For FlIP (Fluid Implicit Particle)
    FloatGrid::Ptr mPressure;
    Int32Grid::Ptr mFlags;
};


SmokeSolver::SmokeSolver(float const voxelSize) : mVoxelSize(voxelSize)
{
    initialize();
}


void
SmokeSolver::initialize() {
    using BBox = math::BBox<Vec3s>;

    mXform = math::Transform::createLinearTransform(mVoxelSize);
    float const padding = 2 * mVoxelSize;

    Vec3s minFI = Vec3s(0.f, 0.f, 0.f);
    Vec3s maxFI = Vec3s(2.f + mVoxelSize, 4.f + mVoxelSize, 5.f + mVoxelSize);
    Coord minFICoord = mXform->worldToIndexNodeCentered(minFI);
    Coord maxFICoord = mXform->worldToIndexNodeCentered(maxFI);
    FloatGrid::Ptr fluidLSInit = FloatGrid::create(/*bg = */0.f);
    fluidLSInit->denseFill(CoordBBox(minFICoord, maxFICoord), /*value = */ 1.0, /*active = */ true);
    fluidLSInit->setTransform(mXform);

    Vec3s maxIntr = Vec3s(14.f + mVoxelSize, 5.f + mVoxelSize, 5.f + mVoxelSize);
    Coord maxFIIntrCoord = mXform->worldToIndexNodeCentered(maxIntr);
    FloatGrid::Ptr negativeSpace = FloatGrid::create(/*bg = */0.f);
    negativeSpace->denseFill(CoordBBox(minFICoord, maxFIIntrCoord), /*value = */ 1.0, /*active = */ true);
    negativeSpace->setTransform(mXform);

    Vec3s minBBoxvec = Vec3s(-padding, -padding, -padding);
    Vec3s maxBBoxvec = Vec3s(14.f + padding + mVoxelSize, 5.f + padding + mVoxelSize, 5.f + padding + mVoxelSize);
    Coord minBBoxcoord = mXform->worldToIndexNodeCentered(minBBoxvec);
    Coord maxBBoxcoord = mXform->worldToIndexNodeCentered(maxBBoxvec);
    mBBoxLS = FloatGrid::create(/*bg = */0.f);
    mBBoxLS->denseFill(CoordBBox(minBBoxcoord, maxBBoxcoord), /*value = */ 1.0, /*active = */ true);
    mBBoxLS->setTransform(mXform);
    mBBoxLS->topologyDifference(*negativeSpace);
    mBBoxLS->topologyDifference(*fluidLSInit);
    mBBoxLS->setName("collider");
    openvdb::tools::pruneInactive(mBBoxLS->tree());

    mPoints = points::denseUniformPointScatter(*fluidLSInit, mPointsPerVoxel);
    mPoints->setName("Points");
    points::appendAttribute<Vec3s>(mPoints->tree(),
                                   "velocity" /* attribute name */,
                                   Vec3s(0.f, 0.f, 0.f) /* uniform value */,
                                   1 /* stride or total count */,
                                   true /* constant stride */,
                                   nullptr /* default value */,
                                   false /* hidden */,
                                   false /* transient */);
    points::appendAttribute<Vec3s>(mPoints->tree(),
                                   "v_pic" /* attribute name */,
                                   Vec3s(0.f, 0.f, 0.f) /* uniform value */,
                                   1 /* stride or total count */,
                                   true /* constant stride */,
                                   nullptr /* default value */,
                                   false /* hidden */,
                                   false /* transient */);
    points::appendAttribute<Vec3s>(mPoints->tree(),
                                   "v_flip" /* attribute name */,
                                   Vec3s(0.f, 0.f, 0.f) /* uniform value */,
                                   1 /* stride or total count */,
                                   true /* constant stride */,
                                   nullptr /* default value */,
                                   false /* hidden */,
                                   false /* transient */);

    openvdb::Index64 count = openvdb::points::pointCount(mPoints->tree());
    std::cout << "PointCount=" << count << std::endl;
}


void
SmokeSolver::particlesToGrid(){
    TreeBase::Ptr baseVTree = points::rasterizeTrilinear<true /* staggered */, Vec3s>(mPoints->tree(), "velocity");

    Vec3STree::Ptr velTree = DynamicPtrCast<Vec3STree>(baseVTree);
    mVCurr = Vec3SGrid::create(velTree);
    mVCurr->setGridClass(GRID_STAGGERED);
    mVCurr->setTransform(mXform);
    mVCurr->setName("v_curr");

    mVNext = Vec3SGrid::create(Vec3s(0.f, 0.f, 0.f));
    (mVNext->tree()).topologyUnion(mVCurr->tree());
    mVNext->setGridClass(GRID_STAGGERED);
    mVNext->setTransform(mXform);
    mVNext->setName("v_next");
}


void
SmokeSolver::addGravity(float const dt) {
    tree::LeafManager<Vec3STree> r(mVCurr->tree());
    SmokeSolver::ApplyGravityOp op(dt, mGravity);
    r.foreach(op);
}


void
SmokeSolver::computeFlipVelocity(float const dt) {
    mVDiff = Vec3SGrid::create(Vec3s(0.f, 0.f, 0.f));
    (mVDiff->tree()).topologyUnion(mVCurr->tree());
    mVDiff->setGridClass(GRID_STAGGERED);
    mVDiff->setTransform(mXform);

    tree::LeafManager<Vec3STree> r(mVDiff->tree());
    SmokeSolver::ComputeFlipVelocityOp op(mVCurr, mVNext, dt, mGravity);
    r.foreach(op);
}


void
SmokeSolver::velocityBCCorrection(Vec3SGrid& vecGrid) {
    auto acc = vecGrid.getAccessor();
    auto bboxAcc = mBBoxLS->getAccessor();

    for (auto iter = vecGrid.beginValueOn(); iter; ++iter) {
        math::Coord ijk = iter.getCoord();
        math::Coord im1jk = ijk.offsetBy(-1, 0, 0);
        math::Coord ip1jk = ijk.offsetBy(1, 0, 0);
        math::Coord ijm1k = ijk.offsetBy(0, -1, 0);
        math::Coord ijp1k = ijk.offsetBy(0, 1, 0);
        math::Coord ijkm1 = ijk.offsetBy(0, 0, -1);
        math::Coord ijkp1 = ijk.offsetBy(0, 0, 1);

        if (bboxAcc.isValueOn(im1jk) || bboxAcc.isValueOn(ip1jk)) {
            auto val = acc.getValue(ijk);
            Vec3s newVal = Vec3s(0, val[1], val[2]);
            acc.setValue(ijk, newVal);
        }
        if (bboxAcc.isValueOn(ijm1k) || bboxAcc.isValueOn(ijp1k)) {
            auto val = acc.getValue(ijk);
            Vec3s newVal = Vec3s(val[0], 0, val[2]);
            acc.setValue(ijk, newVal);
        }
        if (bboxAcc.isValueOn(ijkm1) || bboxAcc.isValueOn(ijkp1)) {
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
    domainMaskGrid->topologyDifference(*mBBoxLS);

    math::pcg::State state = math::pcg::terminationDefaults<ValueType>();
    state.iterations = 100000;
    state.relativeError = state.absoluteError = epsilon;
    SmokeSolver::BoundaryOp bop(mVoxelSize, mBBoxLS, mVCurr);
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
    computeFlipVelocity(dt);
}


void
SmokeSolver::substep(float const dt) {
    particlesToGrid();
    gridVelocityUpdate(dt);
    gridToParticles();
    updateParticles(dt);
}


void
SmokeSolver::updateParticlesVelocity() {
    // Create a leaf iterator for the PointDataTree.
    auto leafIter = (mPoints->tree()).beginLeaf();

    // Retrieve the index from the descriptor.
    // Used to get the array attribute in the functor.
    auto descriptor = leafIter->attributeSet().descriptor();
    Index64 velIdx = descriptor.find("velocity");
    Index64 vPicIdx = descriptor.find("v_pic");
    Index64 vFlipIdx = descriptor.find("v_flip");

    // PIC/FLIP update
    tree::LeafManager<points::PointDataTree> leafManager(mPoints->tree());
    SmokeSolver::FlipUpdateOp op(velIdx, vPicIdx, vFlipIdx, 0.05 /* alpha in PIC/FlIP update */);
    tbb::parallel_for(leafManager.leafRange(), op);
}


void
SmokeSolver::updateParticles(float const dt) {
    updateParticlesVelocity();
    advectParticles(dt);
}


void
SmokeSolver::render() {
    float const dt = 1.f/24.f;
    for (int frame = 0; frame < 200; ++frame) {
        std::cout << "\nframe = " << frame << "\n";
        substep(dt);
        writeVDBs(frame);
        writeVDBsVerbose(frame);
    }
}


void
SmokeSolver::gridToParticles() {
    // Interpolate PIC velocity
    points::boxSample(*mPoints, *mVNext, "v_pic");

    // Interpolate FLIP velocity
    points::boxSample(*mPoints, *mVDiff, "v_flip");
}


void
SmokeSolver::advectParticles(float const dt) {
    Index const integrationOrder = 1;
    int const steps = 1;

    points::advectPoints(*mPoints, *mVNext, integrationOrder, dt, steps);
}


void
SmokeSolver::writeVDBs(int const frame) {
    std::ostringstream ss;
    ss << "smoke_" << std::setw(3) << std::setfill('0') << frame << ".vdb";
    std::string fileName(ss.str());
    io::File file(fileName.c_str());
    file.write({mPoints});
    file.close();
}


void
SmokeSolver::writeVDBsVerbose(int const frame) {
    std::ostringstream ss;
    ss << "smoke_volume_" << std::setw(3) << std::setfill('0') << frame << ".vdb";
    std::string fileName(ss.str());
    openvdb::io::File file(fileName.c_str());

    openvdb::GridPtrVec grids;
    grids.push_back(mBBoxLS);
    grids.push_back(mCollider);
    grids.push_back(mVCurr);
    grids.push_back(mVNext);
    grids.push_back(mDivBefore);
    grids.push_back(mDivAfter);
    grids.push_back(mPressure);
    // grids.push_back(mInterior);
    file.write(grids);
    file.close();
}
}