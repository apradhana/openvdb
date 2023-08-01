// Copyright Contributors to the OpenVDB Project
// SPDX-License-Identifier: MPL-2.0

#include <algorithm>
#include <fstream>
#include <iostream>
#include <math.h>
#include <cmath>
#include <string>
#include <vector>

#include <openvdb/openvdb.h>
#include <openvdb/points/PointConversion.h>
#include <openvdb/points/PointCount.h>
#include <openvdb/util/logging.h>
#include <openvdb/tools/FastSweeping.h>
#include <openvdb/tools/Composite.h> // for tools::compMax
#include <openvdb/tools/VolumeAdvect.h> // for tools::VolumeAdvection
#include <openvdb/tools/PoissonSolver.h> // for poisson solve
#include <openvdb/tree/NodeManager.h> // for post processing bool grid
#include <openvdb/tree/LeafManager.h> // for LeafManager
#include <openvdb/tools/GridOperators.h> // for divergence and gradient

#include <openvdb/points/PointConversion.h>
#include <openvdb/points/PointAttribute.h> // for appendAttribute
#include <openvdb/points/PointCount.h>
#include <openvdb/points/PointScatter.h>
#include <openvdb/points/PointDataGrid.h>
#include <openvdb/points/PointRasterizeTrilinear.h>
#include <openvdb/tools/Morphology.h> // for erodeActiveValues
#include <openvdb/tools/MeshToVolume.h> // for createLevelSetBox
#include <openvdb/points/PointSample.h> // for PointSample
#include <openvdb/points/PointAdvect.h> // for advectPoints
using namespace openvdb;


class FlipSolver {
public:

    FlipSolver(float const voxelSize);

    void render();

private:

    void initializeFreeFall();
    void initializePool();
    void initializeDamBreak();

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
    void writeVDBsDebug(int const frame);
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

    float mVoxelSize = 0.1f;
    Vec3s mGravity = Vec3s(0.f, -9.8f, 0.f);
    int mPointsPerVoxel = 8;
    math::Transform::Ptr mXform;

    points::PointDataGrid::Ptr mPoints;
    FloatGrid::Ptr mBBoxLS;
    FloatGrid::Ptr mCollider;
    FloatGrid::Ptr mDivBefore;
    FloatGrid::Ptr mDivAfter;
    Vec3SGrid::Ptr mVCurr;
    Vec3SGrid::Ptr mVNext;
    Vec3SGrid::Ptr mVDiff; // For FlIP (Fluid Implicit Particle)
    FloatGrid::Ptr mPressure;
    Int32Grid::Ptr mFlags;
    BoolGrid::Ptr mInterior;
};


FlipSolver::FlipSolver(float const voxelSize) : mVoxelSize(voxelSize)
{
    initializeDamBreak();
}


void
FlipSolver::initializeFreeFall() {
    using BBox = math::BBox<Vec3s>;

    mXform = math::Transform::createLinearTransform(mVoxelSize);
    
    auto wsFluidInit = BBox(Vec3s(3.f, 3.f, 3.f) /* min */, Vec3s(4.f, 4.f, 4.f) /* max */);
    FloatGrid::Ptr fluidLSInit = tools::createLevelSetBox<FloatGrid>(wsFluidInit, *mXform);

    auto wsDomain = BBox(Vec3s(0.f, 0.f, 0.f) /* min */, Vec3s(14.f, 0.5f, 14.f) /* max */); // world space domain
    mBBoxLS = tools::createLevelSetBox<FloatGrid>(wsDomain, *mXform);
    mBBoxLS->setGridClass(GRID_LEVEL_SET);
    mBBoxLS->setName("collider");

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
FlipSolver::initializePool() {
    using BBox = math::BBox<Vec3s>;

    mXform = math::Transform::createLinearTransform(mVoxelSize);

    auto cldrBox = BBox(Vec3s(105.f, 100.f, 10.5f) /* min */, Vec3s(107.f, 105.f, 103.5f) /* max */);
    mCollider = tools::createLevelSetBox<FloatGrid>(cldrBox, *mXform);
    mCollider->setGridClass(GRID_LEVEL_SET);
    mCollider->setName("collider");
    
    // auto wsFluidInit = BBox(Vec3s(0.f, 0.f, 0.f) /* min */, Vec3s(3.f * 0.5f, 4.f * 0.5f, 5.f * 0.5f) /* max */);
    // auto wsFluidInit = BBox(Vec3s(2.f, 2.f, 2.f) /* min */, Vec3s(2.1f, 2.1f, 2.1f) /* max */);

    Vec3s minFI = Vec3s(2.f, 2.f, 2.f);
    Vec3s maxFI = Vec3s(3.f, 2.5f, 3.f);
    Vec3s maxFI2 = Vec3s(3.f, 5.1f, 3.f);
    Coord minFIcoord = mXform->worldToIndexNodeCentered(minFI);
    Coord maxFIcoord = mXform->worldToIndexNodeCentered(maxFI);
    Coord maxFIcoord2 = mXform->worldToIndexNodeCentered(maxFI2);
    Vec3s minBBoxvec = Vec3s(1.8f, 1.8f, 1.8f);
    Vec3s maxBBoxvec = Vec3s(3.2f, 3.2f, 3.2f);
    Coord minBBoxcoord = mXform->worldToIndexNodeCentered(minBBoxvec);
    Coord maxBBoxcoord = mXform->worldToIndexNodeCentered(maxBBoxvec);
    auto wsFluidInit = BBox( minFI/* min */,  maxFI/* max */);
    FloatGrid::Ptr fluidLSInit = FloatGrid::create(/*bg = */0.f);
    fluidLSInit->denseFill(CoordBBox(minFIcoord, maxFIcoord), /*value = */ 1.0, /*active = */ true);
    fluidLSInit->setTransform(mXform);
    FloatGrid::Ptr fluidLSInit2 = FloatGrid::create(/*bg = */0.f);
    fluidLSInit2->denseFill(CoordBBox(minFIcoord, maxFIcoord2), /*value = */ 1.0, /*active = */ true);
    fluidLSInit2->setTransform(mXform);

    mBBoxLS = FloatGrid::create(/*bg = */0.f);
    mBBoxLS->denseFill(CoordBBox(minBBoxcoord, maxBBoxcoord), /*value = */ 1.0, /*active = */ true);
    mBBoxLS->setTransform(mXform);
    mBBoxLS->topologyDifference(*fluidLSInit2);
    mBBoxLS->setName("bbox_ls");
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

    openvdb::Index64 count = openvdb::points::pointCount(mPoints->tree());
    std::cout << "PointCount=" << count << std::endl;
}


void
FlipSolver::initializeDamBreak() {
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
FlipSolver::particlesToGrid(){
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
FlipSolver::addGravity(float const dt) {
    auto vCurrAcc = mVCurr->getAccessor();
    auto bboxAcc = mBBoxLS->getAccessor();

    for (auto iter = mVCurr->beginValueOn(); iter; ++iter) {
        auto ijk = iter.getCoord();
        Vec3s newVel = /*bboxAcc.isValueOn(ijk) ? Vec3s(0, 0, 0) : */ vCurrAcc.getValue(ijk) + dt * mGravity;
        vCurrAcc.setValue(ijk, newVel);
    }
}


void
FlipSolver::computeFlipVelocity(float const dt) {
    mVDiff = Vec3SGrid::create(Vec3s(0.f, 0.f, 0.f));
    (mVDiff->tree()).topologyUnion(mVCurr->tree());
    mVDiff->setGridClass(GRID_STAGGERED);
    mVDiff->setTransform(mXform);

    auto vCurrAcc = mVCurr->getAccessor();
    auto vNextAcc = mVNext->getAccessor();
    auto vDiffAcc = mVDiff->getAccessor();

    for (auto iter = mVNext->beginValueOn(); iter; ++iter) {
        auto ijk = iter.getCoord();
        Vec3s val = vNextAcc.getValue(ijk) - vCurrAcc.getValue(ijk) - dt * mGravity;
        vDiffAcc.setValue(ijk, val);
    }
}


void
FlipSolver::velocityBCCorrection(Vec3SGrid& vecGrid) {
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
FlipSolver::pressureProjection(bool print) {
    std::cout << "pressure projection begins" << std::endl;
    using TreeType = FloatTree;
    using ValueType = TreeType::ValueType;
    using MaskGridType = BoolGrid;
    using PCT = openvdb::math::pcg::JacobiPreconditioner<openvdb::tools::poisson::LaplacianMatrix>;

    const ValueType zero = zeroVal<ValueType>();
    const double epsilon = math::Delta<ValueType>::value();

    BoolTree::Ptr interiorMask(new BoolTree(false));
    interiorMask->topologyUnion(mVCurr->tree());
    tools::erodeActiveValues(*interiorMask, /*iterations=*/1, tools::NN_FACE, tools::IGNORE_TILES);
    BoolGrid::Ptr interiorGrid = BoolGrid::create(interiorMask);
    interiorGrid->setTransform(mXform);
    mInterior = interiorGrid->copy();
    mInterior->setName("interior");

    mDivBefore = tools::divergence(*mVCurr);
    mDivBefore->setName("div_before");
    // NOTE: not doing topology intersection
    // (mDivBefore->tree()).topologyIntersection(interiorGrid->tree());
    float divBefore = 0.f;
    auto divAcc = mDivBefore->getAccessor();
    for (auto iter = mDivBefore->beginValueOn(); iter; ++iter) {
        math::Coord ijk = iter.getCoord();
        auto val = divAcc.getValue(ijk);
        if (std::abs(val) > std::abs(divBefore)) {
            divBefore = val;
        }
    }
    std::cout << "\t== divergence before " << divBefore << std::endl;

    MaskGridType* domainMaskGrid = new MaskGridType(*mDivBefore); // match input grid's topology
    domainMaskGrid->topologyDifference(*mBBoxLS);
    (domainMaskGrid->tree()).topologyIntersection(interiorGrid->tree());

    math::pcg::State state = math::pcg::terminationDefaults<ValueType>();
    state.iterations = 100000;
    state.relativeError = state.absoluteError = epsilon;
    FlipSolver::BoundaryOp bop(mVoxelSize, mBBoxLS, mVCurr);

    util::NullInterrupter interrupter;

    FloatTree::Ptr fluidPressure = tools::poisson::solveWithBoundaryConditionsAndPreconditioner<PCT>(
        mDivBefore->tree(), domainMaskGrid->tree(), bop, state, interrupter, /*staggered=*/true);
    FloatGrid::Ptr fluidPressureGrid = FloatGrid::create(fluidPressure);
    fluidPressureGrid->setTransform(mXform);

    auto pressureAcc = fluidPressureGrid->getAccessor();
    for (auto iter = fluidPressureGrid->beginValueOn(); iter; ++iter) {
        math::Coord ijk = iter.getCoord();

        auto val = pressureAcc.getValue(ijk);
        auto divijk = divAcc.getValue(ijk);
        if (print)
        std::cout << "pressure " << ijk << " = " << val << " div = " << divijk << std::endl;
    }

    // Note: need to dilate in order to do one-sided difference
    // because we use a staggered grid velocity field.
    tools::dilateActiveValues(*fluidPressure, /*iterations=*/1, tools::NN_FACE, tools::IGNORE_TILES);
    mPressure = fluidPressureGrid->copy();
    mPressure->setName("pressure");
    // (fluidPressureGrid->tree()).topologyIntersection(interiorGrid->tree());

    auto vCurrAcc = mVCurr->getAccessor();
    auto vNextAcc = mVNext->getAccessor();
    auto boolAcc = interiorGrid->getAccessor();

    int count = 0;
    for (auto iter = mVCurr->beginValueOn(); iter; ++iter) {
        math::Coord ijk = iter.getCoord();
        Vec3s gradijk;
        gradijk[0] = pressureAcc.getValue(ijk) - pressureAcc.getValue(ijk.offsetBy(-1, 0, 0));
        gradijk[1] = pressureAcc.getValue(ijk) - pressureAcc.getValue(ijk.offsetBy(0, -1, 0));
        gradijk[2] = pressureAcc.getValue(ijk) - pressureAcc.getValue(ijk.offsetBy(0, 0, -1));

            auto val = vCurrAcc.getValue(ijk) - gradijk * mVoxelSize;
            vNextAcc.setValue(ijk, val);
            // This is only multiplied by mVoxelSize because in the computation of
            // gradijk, I don't divide by mVoxelSize.
            if (print)
            std::cout << "gradijk = " << ijk << " = " << gradijk * mVoxelSize << "\tvnext = " << vNextAcc.getValue(ijk) << std::endl;
    }

    mDivAfter = tools::divergence(*mVNext);
    mDivAfter->setName("div_after");
    (mDivAfter->tree()).topologyIntersection(interiorGrid->tree());
    float divAfter = 0.f;
    auto divAfterAcc = mDivAfter->getAccessor();
    for (auto iter = mDivAfter->beginValueOn(); iter; ++iter) {
        math::Coord ijk = iter.getCoord();
        auto val = divAfterAcc.getValue(ijk);
        if (std::abs(val) > std::abs(divAfter)) {
            divAfter = val;
        }
    }
    std::cout << "\t== divergence after pp = " << divAfter << std::endl;

    std::cout << "Success: " << state.success << std::endl;
    std::cout << "Iterations: " << state.iterations << std::endl;
    std::cout << "Relative error: " << state.relativeError << std::endl;
    std::cout << "Absolute error: " << state.absoluteError << std::endl;
    std::cout << "before dilate solution->activeVoxelCount() =  " << fluidPressure->activeVoxelCount() << std::endl;
    std::cout << "pressure projection ends" << std::endl;
}


void
FlipSolver::gridVelocityUpdate(float const dt) {
    addGravity(dt);
    velocityBCCorrection(*mVCurr);
    pressureProjection(false /* print */);
    velocityBCCorrection(*mVNext);
    computeFlipVelocity(dt);
}


void
FlipSolver::substep(float const dt) {
    particlesToGrid();
    gridVelocityUpdate(dt);
    gridToParticles();
    updateParticles(dt);
}


void
FlipSolver::updateParticlesVelocity() {
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
    FlipSolver::FlipUpdateOp op(velIdx, vPicIdx, vFlipIdx, 0.05 /* alpha in PIC/FlIP update */);
    tbb::parallel_for(leafManager.leafRange(), op);
}


void
FlipSolver::updateParticles(float const dt) {
    updateParticlesVelocity();
    advectParticles(dt);
}


void
FlipSolver::render() {
    float const dt = 1.f/24.f;
    for (int frame = 0; frame < 200; ++frame) {
        std::cout << "frame = " << frame << "\n";
        substep(dt);
        writeVDBs(frame);
        writeVDBsDebug(frame);
    }
}


void
FlipSolver::gridToParticles() {
    // Interpolate PIC velocity
    points::boxSample(*mPoints, *mVNext, "v_pic");

    // Interpolate FLIP velocity
    points::boxSample(*mPoints, *mVDiff, "v_flip");
}


void
FlipSolver::advectParticles(float const dt) {
    Index const integrationOrder = 1;
    int const steps = 1;

    points::advectPoints(*mPoints, *mVNext, integrationOrder, dt, steps);
}


void
FlipSolver::writeVDBs(int const frame) {
    std::ostringstream ss;
    ss << "water_" << std::setw(3) << std::setfill('0') << frame << ".vdb";
    std::string fileName(ss.str());
    io::File file(fileName.c_str());
    file.write({mPoints});
    file.close();
}


void
FlipSolver::writeVDBsDebug(int const frame) {
    std::ostringstream ss;
    ss << "water_volume_" << std::setw(3) << std::setfill('0') << frame << ".vdb";
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
    grids.push_back(mInterior);
    file.write(grids);
    file.close();
}


class SmokeSolver {
public:
    SmokeSolver();

    void initialize();
    void foobar();
    void foobar2();

    void substep(float const dt);

    void render();

    void updateEmitters();

    void pressureProjection();

    void advectVelocity(float const dt);

    void advectDensity(float const dt);

    void vorticityConfinement();

    void writeVDBs(int const frame);

    struct BoundaryFoobarOp {
        void operator()(const openvdb::Coord& ijk, const openvdb::Coord& neighbor,
            double& source, double& diagonal) const
        {
            //if (neighbor.x() == ijk.x() && neighbor.z() == ijk.z()) {
            //    // Workaround for spurious GCC 4.8 -Wstrict-overflow warning:
            //    const openvdb::Coord::ValueType dy = (ijk.y() - neighbor.y());
            //    if (dy > 0) source -= 1.0;
            //    else diagonal -= 1.0;
            //}
        }
    };

private:
    openvdb::FloatGrid::Ptr mDensity;
    std::vector<openvdb::FloatGrid::Ptr> mColliders;
    openvdb::Vec3SGrid::Ptr mVCurr;
    openvdb::Vec3SGrid::Ptr mVNext;
};

SmokeSolver::SmokeSolver() {
    this->initialize();
}

struct BoundaryOp {
    void operator()(const openvdb::Coord& ijk, const openvdb::Coord& neighbor,
        double& source, double& diagonal) const
    {
        // Boundary conditions:
        // (-X) - Dirichlet, sin(y/5)
        // (+X) - Dirichlet, -sin(y/5)
        // (-Y, +Y, -Z, +Z) - Neumann, dp/d* = 0
        //
        // There's nothing to do for zero Neumann
        //
        // This is the -X face of the domain:
        if (neighbor.x() + 1 == ijk.x()) {
            const double bc = sin(ijk.y() * 0.2);
            source -= bc;
            diagonal -= 1.0; // must "add back" the diagonal's contribution for Dirichlet BCs!!!
        }

        // This is the +X face of the domain:
        if (neighbor.x() - 1 == ijk.x()) {
            const double bc = -sin(ijk.y() * 0.2);

            source -= bc;
            diagonal -= 1.0; // must "add back" the diagonal's contribution for Dirichlet BCs!!!
        }
    }
};

void
testPoissonSolve() {
    using TreeType = FloatTree;
    using ValueType = TreeType::ValueType;

    const int N = 62;
    const ValueType zero = zeroVal<ValueType>();
    const double epsilon = math::Delta<ValueType>::value();
    TreeType source(/*background=*/zero);
    source.fill(CoordBBox(Coord(-N, -N/2, -N/2), Coord(N, N/2-1, N/2-1)), /*value=*/zero);
    math::pcg::State state = math::pcg::terminationDefaults<ValueType>();
    state.iterations = 100;
    state.relativeError = state.absoluteError = epsilon;
    util::NullInterrupter interrupter;
    typename TreeType::Ptr solution = tools::poisson::solveWithBoundaryConditions(
        source, BoundaryOp(), state, interrupter, /*staggered=*/true);
    std::cout << "Success: " << state.success << std::endl;
    std::cout << "Iterations: " << state.iterations << std::endl;
    std::cout << "Relative error: " << state.relativeError << std::endl;
    std::cout << "Absolute error: " << state.absoluteError << std::endl;
    auto grid = Grid<TreeType>::create(solution);
    grid->setTransform(math::Transform::createLinearTransform(10.0/N));
    io::File file("poisson.vdb");
    file.write({grid});
    file.close();
}

void
testDivergence()
{
    // This test is slightly different than the one above for sanity
    // checking purposes.
    using namespace openvdb;

    Vec3SGrid::Ptr inGrid = VectorGrid::create();
    inGrid->setGridClass(GRID_STAGGERED);
    auto& inTree = inGrid->tree();

    int const GRID_DIM = 10;
    int dim = GRID_DIM;
    for (int x = -dim; x<dim; ++x) {
        for (int y = -dim; y<dim; ++y) {
            for (int z = -dim; z<dim; ++z) {
                inTree.setValue(Coord(x,y,z),
                    VectorTree::ValueType(float(x), float(y), float(z)));
            }
        }
    }

    std::cout << "math::Pow3(2 *dim) = " << math::Pow3(2*dim) << "\t" << int(inTree.activeVoxelCount()) << std::endl;

    FloatGrid::Ptr divGrid = tools::divergence(*inGrid);

    FloatGrid::ConstAccessor accessor = divGrid->getConstAccessor();
    --dim;
    for (int x = -dim; x<dim; ++x) {
        for (int y = -dim; y<dim; ++y) {
            for (int z = -dim; z<dim; ++z) {
                Coord xyz(x,y,z);
                Vec3STree::ValueType v = inTree.getValue(xyz);

                const float d = accessor.getValue(xyz);
                if (std::abs(d-3.0f) > 1.e-5) {
                    std::cout << "error" << std::endl;
                }
            }
        }
    }
}


void
SmokeSolver::initialize() {
    openvdb::io::File fileSrc("/home/andre/dev/openvdb_aswf/_data/sphere_fog.vdb");
    fileSrc.open();
    openvdb::GridBase::Ptr baseGrid;
    auto nameIter = fileSrc.beginName();
    baseGrid = fileSrc.readGrid(nameIter.gridName());
    fileSrc.close();
    mDensity = openvdb::gridPtrCast<openvdb::FloatGrid>(baseGrid);
    mDensity->setName("density");

    openvdb::io::File fileBB("/home/andre/dev/openvdb_aswf/_data/bounding_box.vdb");
    fileBB.open();
    openvdb::GridBase::Ptr bbBaseGrid;
    auto bbNameIter = fileBB.beginName();
    bbBaseGrid = fileBB.readGrid(bbNameIter.gridName());
    fileBB.close();
    openvdb::FloatGrid::Ptr fogBB = openvdb::gridPtrCast<openvdb::FloatGrid>(bbBaseGrid);
    // TODO: reconsider the background of the velocity field
    mVCurr = openvdb::Vec3SGrid::create(openvdb::Vec3s(0.f, 1.f, 0.f));
    mVCurr->setTransform(mDensity->transform().copy());
    mVCurr->setName("velocity");
    mVCurr->setGridClass(GRID_STAGGERED);
    mVCurr->tree().topologyUnion(fogBB->tree());
    auto acc = mVCurr->getAccessor();
    for (openvdb::Vec3SGrid::ValueOnIter iter = mVCurr->beginValueOn(); iter; ++iter) {
        acc.setValue(iter.getCoord(), openvdb::Vec3s(0, 1, 0));
    }

    FloatGrid::Ptr divGrid = tools::divergence(*mVCurr);

    FloatGrid::ConstAccessor accessor = divGrid->getConstAccessor();
    float divSum = 0.f;
    for (openvdb::Vec3SGrid::ValueOnIter iter = mVCurr->beginValueOn(); iter; ++iter) {
        divSum += accessor.getValue(iter.getCoord());
    }
    std::cout << "divSum = " << divSum << std::endl;

}

void
SmokeSolver::updateEmitters() {
    openvdb::io::File fileSrc("/home/andre/dev/openvdb_aswf/_data/sphere_fog.vdb");
    fileSrc.open();
    openvdb::GridBase::Ptr baseGrid;
    auto nameIter = fileSrc.beginName();
    baseGrid = fileSrc.readGrid(nameIter.gridName());
    fileSrc.close();
    openvdb::FloatGrid::Ptr densitySrc = openvdb::gridPtrCast<openvdb::FloatGrid>(baseGrid);
    openvdb::tools::compMax(*mDensity /* result grid */, *densitySrc);
}

void
SmokeSolver::pressureProjection() {
}

void
SmokeSolver::advectDensity(float const dt) {
    using AdvT = openvdb::tools::VolumeAdvection<Vec3fGrid>;
    using SamplerT = openvdb::tools::Sampler<1>;

    AdvT advection(*mVCurr);
    advection.setIntegrator(tools::Scheme::MAC);

    auto const& xform = mDensity->transform();
    auto const minCoord = xform.worldToIndexCellCentered(Vec3f(-.5f, 0.f, -.5f));
    auto const maxCoord = xform.worldToIndexCellCentered(Vec3f(.5f, 3.f, .5f));
    BoolGrid::Ptr mask = BoolGrid::create(false);
    mask->fill(CoordBBox(minCoord, maxCoord), true);
    mask->setTransform(mDensity->transform().copy());

    auto newDensity = advection.advect<FloatGrid, BoolGrid, SamplerT>(*mDensity, *mask, dt);
    mDensity = newDensity;
}

void
SmokeSolver::advectVelocity(float const dt) {

}

void
SmokeSolver::vorticityConfinement() {

}

void
SmokeSolver::substep(float const dt) {
    this->updateEmitters();
    this->pressureProjection();
    this->advectVelocity(dt);
    this->advectDensity(dt);
    this->vorticityConfinement();
}

void
SmokeSolver::render() {
    float const dt = 0.1f;
    for (int frame = 0; frame < 100; ++frame) {
        std::cout << "frame = " << frame << "\n";
        substep(dt);
        this->writeVDBs(frame);
    }
}

void
SmokeSolver::writeVDBs(int const frame) {
    std::ostringstream ss;
    ss << "smoke_" << std::setw(3) << std::setfill('0') << frame << ".vdb";
    std::string fileName(ss.str());
    openvdb::io::File file(fileName.c_str());

    openvdb::GridPtrVec grids;
    grids.push_back(mDensity);
    grids.push_back(mVCurr);
    file.write(grids);
    file.close();
}


// https://github.com/AcademySoftwareFoundation/openvdb/blob/master/openvdb/openvdb/tools/FastSweeping.h
struct SetTrueNearNarrowBandKernel {
    SetTrueNearNarrowBandKernel(float epsilon, openvdb::FloatGrid::Ptr sdf) : mEpsilon(epsilon), mSdf(sdf) {}

    // Root node, internal nodes, and leaf nodes
    template<typename NodeT>
    void operator()(NodeT& node, size_t = 1) const
    {
        openvdb::FloatGrid::Accessor acc = mSdf->getAccessor();
        for (auto iter = node.beginValueAll(); iter; ++iter) {
            auto const ijk = iter.getCoord();
            if (std::abs(acc.getValue(ijk)) < mEpsilon) {
                iter.setValue(true);
            } else {
                iter.setValue(false);
            }
        }
    }

  float mEpsilon;
  openvdb::FloatGrid::Ptr mSdf;
};// SetTrueNearNarrowBandKernel


struct SetZeroNearNarrowBandKernel {
    SetZeroNearNarrowBandKernel(float epsilon, openvdb::FloatGrid::Ptr sdf) : mEpsilon(epsilon), mSdf(sdf) {}

    // Root node, internal nodes, and leaf nodes
    template<typename NodeT>
    void operator()(NodeT& node, size_t = 1) const
    {
        openvdb::FloatGrid::Accessor acc = mSdf->getAccessor();
        for (auto iter = node.beginValueAll(); iter; ++iter) {
            auto const ijk = iter.getCoord();
            if (std::abs(acc.getValue(ijk)) < mEpsilon) {
                iter.setValue(0.0f);
            } else {
                iter.setValue(1.0);
            }
        }
    }

  float mEpsilon;
  openvdb::FloatGrid::Ptr mSdf;
};// SetZeroNearNarrowBandKernel

void convertToBool() {
    openvdb::io::File fileSrc("/home/andre/dev/openvdb_aswf/_data/bunny.vdb");
    fileSrc.open();
    openvdb::GridBase::Ptr baseGrid;
    openvdb::io::File::NameIterator nameIter = fileSrc.beginName();
    baseGrid = fileSrc.readGrid(nameIter.gridName());
    fileSrc.close();
    openvdb::FloatGrid::Ptr sdf = openvdb::gridPtrCast<openvdb::FloatGrid>(baseGrid);

    auto voxelSize = sdf->voxelSize();

    // Create a boolgrid that is true near the zero iso-contour of the level-set
    openvdb::BoolGrid::Ptr bGrid = openvdb::BoolGrid::create(false);
    bGrid->tree().topologyUnion(sdf->tree());
    openvdb::tree::NodeManager<openvdb::BoolTree> nodeManager(bGrid->tree());
    SetTrueNearNarrowBandKernel op(voxelSize.length(), sdf);
    nodeManager.foreachTopDown(op, true /* = threaded*/, 1 /* = grainSize*/);
    bGrid->setTransform(sdf->transform().copy());

    // Write boolgrid to a file
    openvdb::io::File file("boolgrid.vdb");
    openvdb::GridPtrVec grids;
    grids.push_back(bGrid);
    file.write(grids);
    file.close();
}


void convertToOnesAndZeros() {
    openvdb::io::File fileSrc("/home/andre/dev/openvdb_aswf/_data/bunny.vdb");
    fileSrc.open();
    openvdb::GridBase::Ptr baseGrid;
    openvdb::io::File::NameIterator nameIter = fileSrc.beginName();
    baseGrid = fileSrc.readGrid(nameIter.gridName());
    fileSrc.close();
    openvdb::FloatGrid::Ptr sdf = openvdb::gridPtrCast<openvdb::FloatGrid>(baseGrid);

    auto voxelSize = sdf->voxelSize();

    // Create a boolgrid that is true near the zero iso-contour of the level-set
    openvdb::FloatGrid::Ptr grid = openvdb::FloatGrid::create(1.0f);
    grid->tree().topologyUnion(sdf->tree());
    openvdb::tree::NodeManager<openvdb::FloatTree> nodeManager(grid->tree());
    SetZeroNearNarrowBandKernel op(voxelSize.length(), sdf);
    nodeManager.foreachTopDown(op, true /* = threaded*/, 1 /* = grainSize*/);
    grid->setTransform(sdf->transform().copy());

    // Write boolgrid to a file
    openvdb::io::File file("zero_one_grid.vdb");
    openvdb::GridPtrVec grids;
    grids.push_back(grid);
    file.write(grids);
    file.close();
}

void simpleGettingNodes() {
    using TreeType = typename openvdb::FloatGrid::TreeType;
    using RootNodeType = typename TreeType::RootNodeType;
    using NodeChainType = typename RootNodeType::NodeChainType;
    using InternalUpperNodeType = typename NodeChainType::template Get<2>; // Upper internal node
    using InternalLowerNodeType = typename NodeChainType::template Get<1>; // Lower internal node
    using LeafNodeType = typename NodeChainType::template Get<0>; // Leaf node

    openvdb::FloatGrid::Ptr grid = openvdb::FloatGrid::create(10.0);
    openvdb::FloatGrid::Accessor accessor = grid->getAccessor();
    openvdb::Coord xyz(0, 0, 0);
    accessor.setValue(xyz, 1.0);

    std::vector<const InternalLowerNodeType *> lowerNodes;
    std::vector<const InternalUpperNodeType *> upperNodes;
    std::vector<const LeafNodeType *> leafNodes;
    (grid->tree()).getNodes(upperNodes);
    (grid->tree()).getNodes(lowerNodes);
    (grid->tree()).getNodes(leafNodes);
    std::cout << "upperNodes.size() = " << upperNodes.size() << "\tupperNodes[0]->LOG2DIM = " << upperNodes[0]->LOG2DIM << std::endl;
    std::cout << "lowerNodes.size() = " << lowerNodes.size() << "\tlowerNodes[0]->LOG2DIM = " << lowerNodes[0]->LOG2DIM << std::endl;
    std::cout << "leafNodes.size() = " << leafNodes.size() << "\tleafNodes[0]->LOG2DIM = " << leafNodes[0]->LOG2DIM << std::endl;
}

void createUnitBox() {
    using std::sin;
    float const voxelSize = 0.05;
    openvdb::math::Transform::Ptr xform = openvdb::math::Transform::createLinearTransform(voxelSize);
    Vec3f const minWorld(0.0, 0.0, 0.0);
    Vec3f const maxWorld(1.0, 1.0, 1.0);
    Coord minIndex = xform->worldToIndexCellCentered(minWorld);
    Coord maxIndex = xform->worldToIndexCellCentered(maxWorld);

    FloatGrid::Ptr floatGrid = FloatGrid::create(/*bg = */0.0);
    floatGrid->denseFill(CoordBBox(minIndex, maxIndex), /*value = */ 1.0, /*active = */ true);
    floatGrid->setTransform(xform->copy());
    floatGrid->setName("grid");

    auto acc = floatGrid->getAccessor();
    for (auto iter = floatGrid->beginValueOn(); iter; ++iter) {
        auto ijk = iter.getCoord();
        auto p = xform->indexToWorld(ijk);
        auto x = p[0]; auto y = p[1]; auto z = p[2];

        std::cout << "worldspace = " << xform->indexToWorld(ijk) << std::endl;
    }

    // Write boolgrid to a file
    openvdb::io::File file("floatgrid.vdb");
    openvdb::GridPtrVec grids;
    grids.push_back(floatGrid);
    file.write(grids);
    file.close();
}

// Wrapper class for creating a all the VDBs needed for the input of
// using the Heat Method.
struct HeatMethod {
    using Coord = openvdb::Coord;
    using Int32Grid = openvdb::Int32Grid;
    using Transform = openvdb::math::Transform;
    using File = openvdb::io::File;
    using GridPtrVec = openvdb::GridPtrVec;

    HeatMethod(float voxelSize) :
        mVoxelSize(voxelSize)
    {}

    // Returns an int grid for the flags.
    // Also write the grid to a file.
    // Dirichlet = 4
    // Interior = 1
    // Neumann = 0
    openvdb::Int32Grid::Ptr createFlags()
    {
        Transform::Ptr xform = openvdb::math::Transform::createLinearTransform(mVoxelSize);
        Coord minIndex = xform->worldToIndexCellCentered(mMinBBox);
        Coord maxIndex = xform->worldToIndexCellCentered(mMaxBBox);

        Int32Grid::Ptr grid = Int32Grid::create(/*bg = */0.0);
        grid->denseFill(CoordBBox(minIndex, maxIndex), /*value = interior = */ 1, /*active = */ true);
        for (auto iter = grid->beginValueOn(); iter; ++iter) {
            auto ijk = iter.getCoord();
            auto p = xform->indexToWorld(ijk);
            auto x = p[0]; auto y = p[1]; auto z = p[2];

            if (std::abs(x) < mEps || std::abs(1 - x) < mEps
                || std::abs(y) < mEps || std::abs(1 - y) < mEps
                || std::abs(z) < mEps || std::abs(1 - z) < mEps) {
                    iter.setValue(4);
            }
        }

        grid->setTransform(xform->copy());
        grid->setName("flags");

        // Write boolgrid to a file
        openvdb::io::File file("flagsHeat.vdb");
        openvdb::GridPtrVec grids;
        grids.push_back(grid);
        file.write(grids);
        file.close();

        return grid;
    }

    // Returns a float grid with the boundary condition.
    openvdb::FloatGrid::Ptr createDirichletBC()
    {
        using std::sin;

        openvdb::math::Transform::Ptr xform = openvdb::math::Transform::createLinearTransform(mVoxelSize);
        Coord minIndex = xform->worldToIndexCellCentered(mMinBBox);
        Coord maxIndex = xform->worldToIndexCellCentered(mMaxBBox);

        FloatGrid::Ptr grid = FloatGrid::create(/*bg = */0.0);
        grid->denseFill(CoordBBox(minIndex, maxIndex), /*value = */ 0, /*active = */ true);
        grid->setTransform(xform->copy());
        grid->setName("dirichlet_bc");

        for (auto iter = grid->beginValueOn(); iter; ++iter) {
            auto ijk = iter.getCoord();
            auto p = xform->indexToWorld(ijk);
            auto x = p[0]; auto y = p[1]; auto z = p[2];

            if (std::abs(x) < mEps || std::abs(1 - x) < mEps
                || std::abs(y) < mEps || std::abs(1 - y) < mEps
                || std::abs(z) < mEps || std::abs(1 - z) < mEps) {
                    float const val = sin(M_PI * (x + y + z) / 3.0) + x * y * z;
                    iter.setValue(val);
            }
        }
        // Write boolgrid to a file
        openvdb::io::File file("dirichletBCHeat.vdb");
        openvdb::GridPtrVec grids;
        grids.push_back(grid);
        file.write(grids);
        file.close();

        return grid;
    }

    float mVoxelSize = 0.05;
    float mEps = 1.0e-5;
    Vec3f mMinBBox = Vec3f(0.0, 0.0, 0.0);
    Vec3f mMaxBBox = Vec3f(1.0, 1.0, 1.0);
}; // HeatMethod


void heatmethod() {
    std::cout << "heatmethod begins" << std::endl;
    HeatMethod hm(0.1f);
    hm.createFlags();
    hm.createDirichletBC();
    std::cout << "heatmethod ends" << std::endl;
}

void testFogToSdf() {
    // openvdb::io::File fileSrc("/home/andre/dev/openvdb_aswf/_data/bunny_cloud.vdb");
    openvdb::io::File fileSrc("/home/andre/dev/openvdb_aswf/_build_google_groups/openvdb_examples/vdb_example_fluids/fog_test.vdb");
    fileSrc.open();
    openvdb::GridBase::Ptr baseGrid;
    openvdb::io::File::NameIterator nameIter = fileSrc.beginName();
    baseGrid = fileSrc.readGrid(nameIter.gridName());
    fileSrc.close();
    openvdb::FloatGrid::Ptr fog = openvdb::gridPtrCast<openvdb::FloatGrid>(baseGrid);

    auto voxelSize = fog->voxelSize();
    auto sdf = openvdb::tools::fogToSdf(*fog, 0.5 /* = isoValue */, 500 /* = iter */);
    sdf->setTransform(fog->transform().copy());
    sdf->setName("sdf");

    openvdb::io::File file("bunny_sdf_3.vdb");
    openvdb::GridPtrVec grids;
    grids.push_back(sdf);
    file.write(grids);
    file.close();
}

void
SmokeSolver::foobar() {
    using PCT = math::pcg::JacobiPreconditioner<tools::poisson::LaplacianMatrix>;

    // Density
    openvdb::io::File fileSrc("/home/andre/dev/openvdb_aswf/_data/sphere_fog.vdb");
    fileSrc.open();
    openvdb::GridBase::Ptr baseGrid;
    auto nameIter = fileSrc.beginName();
    baseGrid = fileSrc.readGrid(nameIter.gridName());
    fileSrc.close();
    mDensity = openvdb::gridPtrCast<openvdb::FloatGrid>(baseGrid);
    mDensity->setName("density");
    auto const& xform = mDensity->transform();

    // Velocity field
    openvdb::io::File fileBB("/home/andre/dev/openvdb_aswf/_data/bounding_box.vdb");
    fileBB.open();
    openvdb::GridBase::Ptr bbBaseGrid;
    auto bbNameIter = fileBB.beginName();
    bbBaseGrid = fileBB.readGrid(bbNameIter.gridName());
    fileBB.close();
    openvdb::FloatGrid::Ptr fogBB = openvdb::gridPtrCast<openvdb::FloatGrid>(bbBaseGrid);
    mVCurr = openvdb::Vec3SGrid::create(openvdb::Vec3s(0.f, 0.f, 0.f));
    mVCurr->setName("velocity");
    mVCurr->setGridClass(GRID_STAGGERED);
    mVCurr->setTransform(mDensity->transform().copy());
    mVCurr->tree().topologyUnion(fogBB->tree());
    auto vCurrAcc = mVCurr->getAccessor();
    for (openvdb::Vec3SGrid::ValueOnIter iter = mVCurr->beginValueOn(); iter; ++iter) {
        math::Coord ijk = iter.getCoord();
        math::Vec3s xyz = xform.indexToWorld(ijk);
        vCurrAcc.setValue(ijk, openvdb::Vec3s(xyz[0] * xyz[0]));
    }

    // Divergence
    FloatGrid::Ptr divGrid = tools::divergence(*mVCurr);
    std::cout << "divGrid->voxelSize()[0] = " << divGrid->voxelSize()[0] << std::endl;
    divGrid->setName("divergence");
    FloatGrid::ConstAccessor divAcc = divGrid->getConstAccessor();
    float divSum = 0.f;
    for (openvdb::Vec3SGrid::ValueOnIter iter = mVCurr->beginValueOn(); iter; ++iter) {
        divSum += divAcc.getValue(iter.getCoord());
    }
    std::cout << "divSum = " << divSum << std::endl;

    // Conjugate Gradient
    std::cout << std::endl << " ==== Conjugate Gradient ====" << std::endl;
    math::pcg::State state = math::pcg::terminationDefaults<float>();
    state.iterations = 100;
    state.relativeError = state.absoluteError = math::Delta<float>::value();
    util::NullInterrupter interrupter;
    openvdb::FloatTree::Ptr pressureTree = tools::poisson::solveWithBoundaryConditionsAndPreconditioner<PCT>(
        divGrid->tree(), divGrid->tree(), BoundaryFoobarOp(), state,
        interrupter, /*staggered=*/true);
    openvdb::FloatGrid::Ptr pressureGrid = FloatGrid::create(pressureTree);
    pressureGrid->setTransform(mVCurr->transform().copy());
    pressureGrid->setName("pressure");
    std::cout << "Success: " << state.success << std::endl;
    std::cout << "Iterations: " << state.iterations << std::endl;
    std::cout << "Relative error: " << state.relativeError << std::endl;
    std::cout << "Absolute error: " << state.absoluteError << std::endl;

    // Update vCurr
    std::cout << std::endl << " ==== Pressure Projection ====" << std::endl;
    auto const voxelSize = mVCurr->voxelSize()[0];
    std::cout << "voxelSize = " << voxelSize << std::endl;
    tools::Gradient<FloatGrid> gradientOp(*pressureGrid);
    Vec3SGrid::Ptr gradientOfPressure = gradientOp.process();
    auto pGradAcc = pressureGrid->getAccessor();
    for (openvdb::Vec3SGrid::ValueOnIter iter = mVCurr->beginValueOn(); iter; ++iter) {
        math::Coord ijk = iter.getCoord();
        Vec3s newVel = *iter - voxelSize * voxelSize * pGradAcc.getValue(ijk);
        vCurrAcc.setValue(ijk, newVel);
    }

    // Divergence2
    FloatGrid::Ptr divAfterGrid = tools::divergence(*mVCurr);
    FloatGrid::ConstAccessor divAfterAcc = divAfterGrid->getConstAccessor();
    divSum = 0.f;
    for (openvdb::Vec3SGrid::ValueOnIter iter = mVCurr->beginValueOn(); iter; ++iter) {
        divSum += divAcc.getValue(iter.getCoord());
    }
    std::cout << "divSum after = " << divSum << std::endl;

    // Print:
    // mVCurr
    // divGrid
    // pressureGrid
    openvdb::io::File file("CG_debug.vdb");
    openvdb::GridPtrVec grids;
    grids.push_back(mVCurr);
    grids.push_back(divGrid);
    grids.push_back(pressureGrid);
    file.write(grids);
    file.close();
}


void
SmokeSolver::foobar2() {
    using PCT = math::pcg::JacobiPreconditioner<tools::poisson::LaplacianMatrix>;

    // Reference grid
    openvdb::io::File fileSrc("/home/andre/dev/openvdb_aswf/_data/sphere_fog.vdb");
    fileSrc.open();
    openvdb::GridBase::Ptr baseGrid;
    auto nameIter = fileSrc.beginName();
    baseGrid = fileSrc.readGrid(nameIter.gridName());
    fileSrc.close();
    openvdb::FloatGrid::Ptr refGrid = gridPtrCast<openvdb::FloatGrid>(baseGrid);
    auto const& xform = refGrid->transform();
    float const voxelSize = refGrid->voxelSize()[0];
    std::cout << "voxelSize = " << voxelSize << std::endl;

    // Density
    std::cout << "setting up density" << std::endl;
    mDensity = openvdb::FloatGrid::create(/*background value=*/0.f);
    mDensity->setName("density");
    mDensity->setTransform(refGrid->transform().copy());
    auto dstAcc = mDensity->getAccessor();
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 2; ++j) {
            for (int k = 0; k < 2; ++k) {
                Coord ijk(i, j, k);
                std::cout << "ijk = " << ijk << std::endl;
                dstAcc.setValue(ijk, 1.0);
            }
        }
    }

    // Velocity field
    std::cout << "setting up velocity" << std::endl;
    mVCurr = openvdb::Vec3SGrid::create(openvdb::Vec3s(0.f, 0.f, 0.f));
    mVCurr->setName("velocity");
    mVCurr->setGridClass(GRID_STAGGERED);
    mVCurr->setTransform(refGrid->transform().copy());
    //mVCurr->tree().topologyUnion(fogBB->tree());
    auto vCurrAcc = mVCurr->getAccessor();
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 2; ++j) {
            for (int k = 0; k < 2; ++k) {
                Coord ijk(i, j, k);
                Vec3s xyz(i * voxelSize, j * voxelSize, k * voxelSize);
                vCurrAcc.setValue(ijk, xyz);
            }
        }
    }
    std::cout << "done with velocity" << std::endl;

    // Divergence
    std::cout << "computing divergence" << std::endl;
    FloatGrid::Ptr divGrid = tools::divergence(*mVCurr);
    divGrid->setName("divergence");
    //divGrid->setTransform(refGrid->transform().copy());
    FloatGrid::ConstAccessor divAcc = divGrid->getConstAccessor();
    float divSum = 0.f;
    for (auto iter = divGrid->beginValueOn(); iter; ++iter) {
        Coord ijk = iter.getCoord();
        float divVal = divAcc.getValue(ijk);
        divSum += divVal;
        std::cout << "div val = " << divVal << ", coord = " << ijk << std::endl;
    }
    std::cout << "divSum = " << divSum << std::endl;
    std::cout << "do i even get here??" << std::endl;

    // // Conjugate Gradient
    // std::cout << std::endl << " ==== Conjugate Gradient ====" << std::endl;
    // math::pcg::State state = math::pcg::terminationDefaults<float>();
    // state.iterations = 100;
    // state.relativeError = state.absoluteError = math::Delta<float>::value();
    // util::NullInterrupter interrupter;
    // openvdb::FloatTree::Ptr pressureTree = tools::poisson::solveWithBoundaryConditionsAndPreconditioner<PCT>(
    //     divGrid->tree(), divGrid->tree(), BoundaryFoobarOp(), state,
    //     interrupter, /*staggered=*/true);
    // openvdb::FloatGrid::Ptr pressureGrid = FloatGrid::create(pressureTree);
    // pressureGrid->setTransform(mVCurr->transform().copy());
    // pressureGrid->setName("pressure");
    // std::cout << "Success: " << state.success << std::endl;
    // std::cout << "Iterations: " << state.iterations << std::endl;
    // std::cout << "Relative error: " << state.relativeError << std::endl;
    // std::cout << "Absolute error: " << state.absoluteError << std::endl;

    // // Update vCurr
    // std::cout << std::endl << " ==== Pressure Projection ====" << std::endl;
    // auto const voxelSize = mVCurr->voxelSize()[0];
    // std::cout << "voxelSize = " << voxelSize << std::endl;
    // tools::Gradient<FloatGrid> gradientOp(*pressureGrid);
    // Vec3SGrid::Ptr gradientOfPressure = gradientOp.process();
    // auto pGradAcc = pressureGrid->getAccessor();
    // for (openvdb::Vec3SGrid::ValueOnIter iter = mVCurr->beginValueOn(); iter; ++iter) {
    //     math::Coord ijk = iter.getCoord();
    //     Vec3s newVel = *iter - voxelSize * voxelSize * pGradAcc.getValue(ijk);
    //     vCurrAcc.setValue(ijk, newVel);
    // }

    // // Divergence2
    // FloatGrid::Ptr divAfterGrid = tools::divergence(*mVCurr);
    // FloatGrid::ConstAccessor divAfterAcc = divAfterGrid->getConstAccessor();
    // divSum = 0.f;
    // for (openvdb::Vec3SGrid::ValueOnIter iter = mVCurr->beginValueOn(); iter; ++iter) {
    //     divSum += divAcc.getValue(iter.getCoord());
    // }
    // std::cout << "divSum after = " << divSum << std::endl;

    // Print:
    // mVCurr
    // divGrid
    // pressureGrid
    openvdb::io::File file("CG_debug.vdb");
    openvdb::GridPtrVec grids;
    grids.push_back(mVCurr);
    grids.push_back(divGrid);
    //grids.push_back(pressureGrid);
    file.write(grids);
    file.close();
}


struct SimpleExampleBoundaryOp {

    SimpleExampleBoundaryOp(const float vs) : voxelSize(vs) {
        xform = math::Transform::createLinearTransform(voxelSize);
    }

    void operator()(const openvdb::Coord& ijk, const openvdb::Coord& neighbor,
        double& source, double& diagonal) const {
        using std::sin;
        auto xyzNgbr = xform->indexToWorld(neighbor);
        //double bc = sin(xyzNgbr[0]);//xyzNgbr[0] * xyzNgbr[0] + xyzNgbr[1] * xyzNgbr[1] + xyzNgbr[2] * xyzNgbr[2];
        // double bc = xyzNgbr[0] * xyzNgbr[0] + xyzNgbr[1] * xyzNgbr[1] + xyzNgbr[2] * xyzNgbr[2];
        double bc = sin(2 * M_PI * xyzNgbr[0]);
        // left x-face
        if (neighbor.x() + 1 == ijk.x() /* left x-face */ ||
            neighbor.x() - 1 == ijk.x() /* right x-face */ ||
            neighbor.y() + 1 == ijk.y() /* bottom y-face */ ||
            neighbor.y() - 1 == ijk.y() /* top y-face */ ||
            neighbor.z() + 1 == ijk.z() /* back z-face */ ||
            neighbor.z() - 1 == ijk.z() /* front z-face */) {
            diagonal -= 1.0;
            source -= bc;
        }
    }

    float voxelSize;
    openvdb::math::Transform::Ptr xform;
};


struct SimpleLinearBoundaryOp {

    SimpleLinearBoundaryOp(const float vs) : voxelSize(vs) {
        xform = math::Transform::createLinearTransform(voxelSize);
    }

    void operator()(const openvdb::Coord& ijk, const openvdb::Coord& neighbor,
        double& source, double& diagonal) const
    {
        using std::sin;
        auto xyzNgbr = xform->indexToWorld(neighbor);
        double bc = sin(xyzNgbr[0]);//xyzNgbr[0] * xyzNgbr[0] + xyzNgbr[1] * xyzNgbr[1] + xyzNgbr[2] * xyzNgbr[2];
        // left x-face
        if (neighbor.x() + 1 == ijk.x() /* left x-face */ ||
            neighbor.x() - 1 == ijk.x() /* right x-face */ ||
            neighbor.y() + 1 == ijk.y() /* bottom y-face */ ||
            neighbor.y() - 1 == ijk.y() /* top y-face */ ||
            neighbor.z() + 1 == ijk.z() /* back z-face */ ||
            neighbor.z() - 1 == ijk.z() /* front z-face */) {
            diagonal -= 1.0;
            source -= bc;
        }
    }

    float voxelSize;
    openvdb::math::Transform::Ptr xform;
};

void
checkPoisson() {
    using TreeType = FloatTree;
    using ValueType = TreeType::ValueType;
    using std::sin;

    const ValueType zero = zeroVal<ValueType>();
    const double epsilon = math::Delta<ValueType>::value();
    
    std::vector<int> voxelsPerDims {4, 8, 16, 32, 64, 128};
    
    for (int i = 0; i < voxelsPerDims.size(); ++i) {
        int const N = voxelsPerDims[i];
        float const voxelSize = 1.0f/static_cast<float>(N);
        auto const xform = math::Transform::createLinearTransform(voxelSize);
        
        FloatTree::Ptr source(new FloatTree(0.f));
        source->denseFill(CoordBBox(Coord(1, 1, 1), Coord(N-1, N-1, N-1)), /* value = */0.f);
        FloatGrid::Ptr sourceGrid = Grid<FloatTree>::create(source);
        sourceGrid->setTransform(xform);
        auto srcAcc = sourceGrid->getAccessor();
        
        FloatTree::Ptr trueSln(new FloatTree(0.f));
        trueSln->denseFill(CoordBBox(Coord(1, 1, 1), Coord(N-1, N-1, N-1)), /* value = */0.f);
        FloatGrid::Ptr slnGrid = Grid<FloatTree>::create(trueSln);
        slnGrid->setTransform(xform);
        auto slnAcc = slnGrid->getAccessor();
        
        for (auto iter = sourceGrid->beginValueOn(); iter; ++iter) {
            auto ijk = iter.getCoord();
            auto xyz = xform->indexToWorld(ijk);
            // float sln = sin(2 * M_PI * xyz[0]) + sin(2 * M_PI * xyz[1]) + sin(2 * M_PI * xyz[2]);
            float sln = sin(2 * M_PI * xyz[0]);
            float rhs = -4.f * M_PI * M_PI * sln * voxelSize * voxelSize;
            // float sln = xyz[0] * xyz[0] + xyz[1] * xyz[1] + xyz[2] * xyz[2];
            // float rhs = 6.f * voxelSize * voxelSize;
            srcAcc.setValue(ijk, rhs);
            slnAcc.setValue(ijk, sln);
            // std::cout << "=== ijk" << ijk << " xyz = " << xyz << " = " << srcAcc.getValue(ijk) << std::endl;
        }
        
        SimpleExampleBoundaryOp bcOp(voxelSize);
        math::pcg::State state = math::pcg::terminationDefaults<ValueType>();
        state.iterations = 1000;
        state.relativeError = state.absoluteError = epsilon * 0.00000001;
        util::NullInterrupter interrupter;
        FloatTree::Ptr fluidPressure = tools::poisson::solveWithBoundaryConditions(
                                                                                   sourceGrid->tree(), bcOp, state, interrupter, /*staggered=*/true);
        
        // std::cout << "Success: " << state.success << std::endl;
        // std::cout << "Iterations: " << state.iterations << std::endl;
        // std::cout << "Relative error: " << state.relativeError << std::endl;
        // std::cout << "Absolute error: " << state.absoluteError << std::endl;
        // std::cout << "before dilate solution->activeVoxelCount() =  " << fluidPressure->activeVoxelCount() << std::endl;
        
        double totalError = 0.0;
        double maxError = 0.0;
        
        FloatGrid::Ptr fluidPressureGrid = FloatGrid::create(fluidPressure);
        auto numAcc = fluidPressureGrid->getAccessor();
        for (auto iter = sourceGrid->beginValueOn(); iter; ++iter) {
            auto ijk = iter.getCoord();
            auto xyz = xform->indexToWorld(ijk);
            float err = slnAcc.getValue(ijk) - numAcc.getValue(ijk);
            float vsMult = numAcc.getValue(ijk) * voxelSize;
            float vsSqrMult = numAcc.getValue(ijk) * voxelSize * voxelSize;
            totalError += err * err;
            if (std::abs(err) > 1.0e-5) {
                if (std::abs(err) > maxError) {
                    maxError = std::abs(err);
                }
                // std::cout << "ijk = " << ijk << " xyz = " << xyz << " true sln = " << slnAcc.getValue(ijk) << " ovdb sln = " << vsSqrMult << " err = " << err << std::endl;
                // std::cout << "ijk = " << ijk << " xyz = " << xyz << " true sln = " << slnAcc.getValue(ijk) << " ovdb sln = " << numAcc.getValue(ijk)  << " err = " << err << std::endl;
            }
        }
        double cvgcTest = totalError * voxelSize * voxelSize * voxelSize;
        cvgcTest = std::sqrt(cvgcTest);
        std::cout << "Voxels per dim = " << N << " iterations = " << state.iterations << " cg relative error = " << state.relativeError << " cg absolute error = " << state.absoluteError << " convergence Test = " << cvgcTest << " maxError = max(true sln - num sln) = " << maxError << std::endl;
    }
}


void
checkPoisson2() {
    using TreeType = FloatTree;
    using ValueType = TreeType::ValueType;
    using std::sin;

    const ValueType zero = zeroVal<ValueType>();
    const double epsilon = math::Delta<ValueType>::value();

    int const N = 50;
    float const voxelSize = 1.0f;
    auto const xform = math::Transform::createLinearTransform(voxelSize);

    FloatTree::Ptr source(new FloatTree(0.f));
    source->denseFill(CoordBBox(Coord(1, 1, 1), Coord(N, N, N)), /* value = */0.f);
    FloatGrid::Ptr sourceGrid = Grid<FloatTree>::create(source);
    sourceGrid->setTransform(xform);
    auto srcAcc = sourceGrid->getAccessor();

    FloatTree::Ptr trueSln(new FloatTree(0.f));
    trueSln->denseFill(CoordBBox(Coord(1, 1, 1), Coord(N, N, N)), /* value = */0.f);
    FloatGrid::Ptr slnGrid = Grid<FloatTree>::create(trueSln);
    slnGrid->setTransform(xform);
    auto slnAcc = slnGrid->getAccessor();

    for (auto iter = sourceGrid->beginValueOn(); iter; ++iter) {
        auto ijk = iter.getCoord();
        auto xyz = xform->indexToWorld(ijk);
        // float sln = sin(2 * M_PI * xyz[0]) + sin(2 * M_PI * xyz[1]) + sin(2 * M_PI * xyz[2]);
        float sln = sin(xyz[0]); //xyz[0] * xyz[0] + xyz[1] * xyz[1] + xyz[2] * xyz[2];
        float rhs = -sin(xyz[0]);
        srcAcc.setValue(ijk, rhs);
        slnAcc.setValue(ijk, sln);
        if (ijk[0] == 12 && ijk[1] == 10 && ijk[2] == 10) {
            std::cout << "=== ijk" << ijk << " xyz = " << xyz << " = " << srcAcc.getValue(ijk) << std::endl;
        }
    }

    SimpleLinearBoundaryOp bcOp(voxelSize);
    math::pcg::State state = math::pcg::terminationDefaults<ValueType>();
    state.iterations = 1000000;
    state.relativeError = state.absoluteError = epsilon * 0.00000001;
    util::NullInterrupter interrupter;
    FloatTree::Ptr fluidPressure = tools::poisson::solveWithBoundaryConditions(
        sourceGrid->tree(), bcOp, state, interrupter, /*staggered=*/true);

    std::cout << "Success: " << state.success << std::endl;
    std::cout << "Iterations: " << state.iterations << std::endl;
    std::cout << "Relative error: " << state.relativeError << std::endl;
    std::cout << "Absolute error: " << state.absoluteError << std::endl;
    std::cout << "before dilate solution->activeVoxelCount() =  " << fluidPressure->activeVoxelCount() << std::endl;
    std::cout << "epsilon = " << epsilon << std::endl;

    double totalError = 0.0;

    FloatGrid::Ptr fluidPressureGrid = FloatGrid::create(fluidPressure);
    auto numAcc = fluidPressureGrid->getAccessor();

    std::cout << "active voxel count: sourceGrid = " << sourceGrid->activeVoxelCount() <<
        " solution = " << slnGrid->activeVoxelCount() <<
        " fluid pressure = " << fluidPressureGrid->activeVoxelCount() << std::endl;
    int count = 0;
    double maxError = 0.0;
    for (auto iter = sourceGrid->beginValueOn(); iter; ++iter) {
        auto ijk = iter.getCoord();
        auto xyz = xform->indexToWorld(ijk);
        float err = slnAcc.getValue(ijk) - numAcc.getValue(ijk);
        float num = numAcc.getValue(ijk);
        totalError += err * err;
        if (std::abs(err) > 1.0e-5) {
            count++;
            //std::cout << "dif count = " << count << " ijk = " << ijk << " xyz = " << xyz << " true sln = " << slnAcc.getValue(ijk) << " ovdb sln = " << num << " err = " << err << std::endl;
            if (std::abs(err) > maxError) {
                maxError = std::abs(err);
            }
            //std::cout << "ijk = " << ijk << " xyz = " << xyz << " true sln = " << slnAcc.getValue(ijk) << " ovdb sln = " << numAcc.getValue(ijk)  << " err = " << err << std::endl;
        }
    }
    double cvgcTest = totalError * voxelSize * voxelSize * voxelSize;
    cvgcTest = std::sqrt(cvgcTest);
    std::cout << "convergence Test = " << cvgcTest << " maxError = " << maxError << std::endl;
}

void
simpleFlip() {
    using BBox = math::BBox<Vec3s>;
    int constexpr pointsPerVoxel = 8;
    float const voxelSize = 0.1f;
    float dt = 1;
    Vec3s const gravity(0.f, -9.8f, 0.f);

    math::Transform::Ptr xform = math::Transform::createLinearTransform(voxelSize);
    
    std::cout << "simple flip begins" << std::endl;
    
    auto wsDomain = BBox(Vec3s(0.f, 0.f, 0.f) /* min */, Vec3s(14.f, 5.f, 5.f) /* max */);
    auto wsFluidInit = BBox(Vec3s(0.f, 0.f, 0.f) /* min */, Vec3s(3.f, 4.f, 5.f) /* max */);
    
    FloatGrid::Ptr fluidLSInit = tools::createLevelSetBox<FloatGrid>(wsFluidInit, *xform);
    
    auto points = points::denseUniformPointScatter(*fluidLSInit, pointsPerVoxel);
    points->setName("Points");
    openvdb::io::File("mypoints_orig.vdb").write({points});
    // Get number of points 
    int pointCount = static_cast<int>(points::pointCount(points->tree()));
    // Populate velocity attribute
    points::appendAttribute<Vec3s>(points->tree(), "velocity",
                                   /*uniformValue*/ Vec3s(0.f, 0.f, 0.f),
                                   /*stride or total count=*/1,
                                   /*constantStride=*/true,
                                   /*defaultValue*/nullptr,
                                   /*hidden=*/false, /*transient=*/false);

    TreeBase::Ptr baseVTree = points::rasterizeTrilinear</*staggered=*/true, Vec3s>(points->tree(), "velocity");


    Vec3STree::Ptr velTree = DynamicPtrCast<Vec3STree>(baseVTree);
    Vec3SGrid::Ptr mVCurr = Vec3SGrid::create(velTree);
    mVCurr->setGridClass(GRID_STAGGERED);
    mVCurr->setTransform(xform);

    Vec3SGrid::Ptr mVNext = Vec3SGrid::create(Vec3s(0.f, 0.f, 0.f));
    (mVNext->tree()).topologyUnion(mVCurr->tree());
    mVNext->setGridClass(GRID_STAGGERED);
    mVNext->setTransform(xform);

    auto vCurrAcc = mVCurr->getAccessor();
    auto vNextAcc = mVNext->getAccessor();

    for (auto iter = mVCurr->beginValueOn(); iter; ++iter) {
        auto ijk = iter.getCoord();
        Vec3s newVel = vCurrAcc.getValue(ijk) + dt * gravity;
        vNextAcc.setValue(ijk, newVel);
    }

    // TODO: interpolate back
    points::AttributeHandle<Vec3s>::Ptr velHandle =
        points::AttributeHandle<Vec3s>::create(
            points->tree().cbeginLeaf()->attributeArray("velocity"));
    std::cout << "velHandle.get() = " << velHandle.get() << std::endl;

    // nearest-neighbour staggered sampling
    points::boxSample(*points, *mVNext, "velocity");
    Index const integrationOrder = 1;
    int const steps = 1;
    double const timeStep = 1.0/24.0;
    points::advectPoints(*points, *mVNext, integrationOrder, timeStep, steps);

    openvdb::io::File("mypoints_next.vdb").write({points});

    std::cout << "fluidLSInit = " << fluidLSInit << std::endl;
    
    std::cout << "simple flip ends" << std::endl;
}

// TO BUILD:
// mkdir build
// cd build
// cmake -DOPENVDB_BUILD_EXAMPLES=ON -DOPENVDB_BUILD_VDB_EXAMPLE_FLUIDS=ON ../
// make -j 8
int
main(int argc, char *argv[])
{
    openvdb::initialize();

    //createUnitBox();
    //foobar();
    // convertToOnesAndZeros();
    // SmokeSolver smokeSim;
    // smokeSim.foobar2();

    FlipSolver flipSim(0.1f /* voxel size */);
    flipSim.render();
    // flipSim.particlesToGrid2();
    // flipSim.pressureProjection2();

    // checkPoisson();
    simpleFlip();


    // solver.render();
    // testPoissonSolve();
    // testDivergence();
}
