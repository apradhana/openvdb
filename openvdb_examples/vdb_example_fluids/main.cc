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

class Vector3 {
public:
    float x, y, z;

    Vector3(float _x, float _y, float _z) : x(_x), y(_y), z(_z) {}
};

namespace helper {
template<typename GridType>
struct ToMaskGrid {
    typedef Grid<typename GridType::TreeType::template ValueConverter<ValueMask>::Type> Type;
};
}
/// @brief Compute the gradient of a scalar grid.
template<
    typename InGridT,
    typename MaskGridType = typename helper::ToMaskGrid<InGridT>::Type,
    typename InterruptT = util::NullInterrupter>
class GradientStaggered
{
public:
    typedef InGridT                                         InGridType;
    typedef typename tools::ScalarToVectorConverter<InGridT>::Type OutGridType;

    GradientStaggered(const InGridT& grid, InterruptT* interrupt = nullptr):
        mInputGrid(grid), mInterrupt(interrupt), mMask(nullptr)
    {
    }

    GradientStaggered(const InGridT& grid, const MaskGridType& mask, InterruptT* interrupt = nullptr):
        mInputGrid(grid), mInterrupt(interrupt), mMask(&mask)
    {
    }

    typename OutGridType::Ptr process(bool threaded = true)
    {
        Functor functor(mInputGrid, mMask, threaded, mInterrupt);
        processTypedMap(mInputGrid.transform(), functor);
        if (functor.mOutputGrid) functor.mOutputGrid->setVectorType(VEC_COVARIANT);
        return functor.mOutputGrid;
    }

protected:
    struct Functor
    {
        Functor(const InGridT& grid, const MaskGridType* mask,
            bool threaded, InterruptT* interrupt):
            mThreaded(threaded), mInputGrid(grid), mInterrupt(interrupt), mMask(mask) {}

        template<typename MapT>
        void operator()(const MapT& map)
        {
            typedef math::Gradient<MapT, math::FD_1ST> OpT;
            tools::gridop::GridOperator<InGridType, MaskGridType, OutGridType, MapT, OpT, InterruptT>
                op(mInputGrid, mMask, map, mInterrupt);
            mOutputGrid = op.process(mThreaded); // cache the result
        }

        const bool                 mThreaded;
        const InGridT&             mInputGrid;
        typename OutGridType::Ptr  mOutputGrid;
        InterruptT*                mInterrupt;
        const MaskGridType*        mMask;
    }; // Private Functor

    const InGridT&       mInputGrid;
    InterruptT*          mInterrupt;
    const MaskGridType*  mMask;
}; // end of Gradient class


class FlipSolver {
public:

    FlipSolver(float const voxelSize);

    void render();

private:

    void initialize();
    void initialize2();
    void initialize3();

    void substep(float const dt);

    // Rasterize particle velocity to the grid
    void particlesToGrid();
    void particlesToGrid2();

    // FLIP update: Interpolate the delta of velocity update (v_np1 - v_n)
    // back to the particle
    void gridToParticles();

    // Update particle position based on velocity on the grid
    void advectParticles(float const dt);

    // Make the velocity on the grid to be divergence free
    void pressureProjection();
    void pressureProjection2();
    void pressureProjection3();
    void pressureProjection4();
    void pressureProjection5();
    void gridVelocityUpdate(float const dt);

    void velocityBCCorrection();

    void addGravity(float const dt);

    void writeVDBs(int const frame);
    void writeVDBsDebug(int const frame);
    struct BoundaryOp {
        BoundaryOp(float const voxelSize,
                   FloatGrid::Ptr bBoxLS,
                   FloatGrid::Ptr collider,
                   Vec3SGrid::Ptr vCurr) :
            voxelSize(voxelSize),
            bBoxLS(bBoxLS),
            collider(collider),
            vCurr(vCurr) {}

        void operator()(const openvdb::Coord& ijk,
                        const openvdb::Coord& neighbor,
                        double& source,
                        double& diagonal) const
        {
            float const dirichletBC = 0.f;
            bool isInsideBBox = bBoxLS->tree().isValueOn(neighbor);
            bool isInsideCollider = collider->tree().isValueOn(neighbor);
            auto vNgbr = vCurr->tree().getValue(neighbor);

            // TODO: Fix this:
            if (/* isInsideCollider ||*/  isInsideBBox) {
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
                source += delta /** 0.5 */ / voxelSize;
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
        FloatGrid::Ptr bBoxLS;
        FloatGrid::Ptr collider;
        Vec3SGrid::Ptr vCurr;
    };

    struct BoundaryFooBarOp {
        BoundaryFooBarOp(float const voxelSize) :
            voxelSize(voxelSize) {

        }

        void operator()(const openvdb::Coord& ijk, const openvdb::Coord& neighbor,
            double& source, double& diagonal) const
        {
            // Boundary conditions:
            // (-X, +X, -Y, -Z, +Z): Neumann dp/dn = 0
            // (+Y): Dirichlet p = 0
            // There is nothing to do for zero value Neumann boundary condition.
            //if (ijk.y()+1 == neighbor.y()) {
            //    source -= 0.0;
                diagonal -= 1.0; // / (voxelSize * voxelSize);
            //}
        }

        float voxelSize;
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
    FloatGrid::Ptr mPressure;
    Int32Grid::Ptr mFlags;
    BoolGrid::Ptr mInterior;
};


FlipSolver::FlipSolver(float const voxelSize) : mVoxelSize(voxelSize)
{
    initialize3();
}


void
FlipSolver::initialize() {
    using BBox = math::BBox<Vec3s>;

    mXform = math::Transform::createLinearTransform(mVoxelSize);

    auto wsDomain = BBox(Vec3s(0.f, 0.f, 0.f) /* min */, Vec3s(14.f, 5.f, 5.f) /* max */); // world space domain
    mBBoxLS = tools::createLevelSetBox<FloatGrid>(wsDomain, *mXform);
    mBBoxLS->setGridClass(GRID_LEVEL_SET);
    mBBoxLS->setName("bbox_ls");

    auto cldrBox = BBox(Vec3s(5.f, 0.f, 1.5f) /* min */, Vec3s(7.f, 5.f, 3.5f) /* max */);
    mCollider = tools::createLevelSetBox<FloatGrid>(cldrBox, *mXform);
    mCollider->setGridClass(GRID_LEVEL_SET);
    mCollider->setName("collider");
    
    // auto wsFluidInit = BBox(Vec3s(0.f, 0.f, 0.f) /* min */, Vec3s(3.f * 0.5f, 4.f * 0.5f, 5.f * 0.5f) /* max */);
    // auto wsFluidInit = BBox(Vec3s(2.f, 2.f, 2.f) /* min */, Vec3s(2.1f, 2.1f, 2.1f) /* max */);
    auto wsFluidInit = BBox(Vec3s(2.f, 2.f, 2.f) /* min */, Vec3s(3.f, 3.f, 3.f) /* max */);
    FloatGrid::Ptr fluidLSInit = tools::createLevelSetBox<FloatGrid>(wsFluidInit, *mXform);

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
FlipSolver::initialize2() {
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
FlipSolver::initialize3() {
    using BBox = math::BBox<Vec3s>;

    mXform = math::Transform::createLinearTransform(mVoxelSize);

    auto cldrBox = BBox(Vec3s(105.f, 100.f, 10.5f) /* min */, Vec3s(107.f, 105.f, 103.5f) /* max */);
    mCollider = tools::createLevelSetBox<FloatGrid>(cldrBox, *mXform);
    mCollider->setGridClass(GRID_LEVEL_SET);
    mCollider->setName("collider");
    
    // auto wsFluidInit = BBox(Vec3s(0.f, 0.f, 0.f) /* min */, Vec3s(3.f * 0.5f, 4.f * 0.5f, 5.f * 0.5f) /* max */);
    // auto wsFluidInit = BBox(Vec3s(2.f, 2.f, 2.f) /* min */, Vec3s(2.1f, 2.1f, 2.1f) /* max */);

    Vec3s minFI = Vec3s(2.f, 2.f, 2.f);
    Vec3s maxFI = Vec3s(2.1f, 2.3f, 2.1f);
    Vec3s maxFI2 = Vec3s(2.1f, 5.1f, 2.1f);
    Coord minFIcoord = mXform->worldToIndexNodeCentered(minFI);
    Coord maxFIcoord = mXform->worldToIndexNodeCentered(maxFI);
    Coord maxFIcoord2 = mXform->worldToIndexNodeCentered(maxFI2);
    Vec3s minBBoxvec = Vec3s(1.8f, 1.8f, 1.8f);
    Vec3s maxBBoxvec = Vec3s(2.3f, 3.2f, 2.3f);
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

    mPoints = points::denseUniformPointScatter(*fluidLSInit, 1 /* mPointsPerVoxel */ );
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
FlipSolver::velocityBCCorrection() {
    auto vCurrAcc = mVCurr->getAccessor();
    auto vNextAcc = mVNext->getAccessor();
    auto bboxAcc = mBBoxLS->getAccessor();

    for (auto iter = mVNext->beginValueOn(); iter; ++iter) {
        math::Coord ijk = iter.getCoord();
        math::Coord im1jk = ijk.offsetBy(-1, 0, 0);
        math::Coord ip1jk = ijk.offsetBy(1, 0, 0);
        math::Coord ijm1k = ijk.offsetBy(0, -1, 0);
        math::Coord ijp1k = ijk.offsetBy(0, 1, 0);
        math::Coord ijkm1 = ijk.offsetBy(0, 0, -1);
        math::Coord ijkp1 = ijk.offsetBy(0, 0, 1);

        if (bboxAcc.isValueOn(im1jk) || bboxAcc.isValueOn(ip1jk)) {
            auto val = vNextAcc.getValue(ijk);
            Vec3s newVal = Vec3s(0, val[1], val[2]);
            vNextAcc.setValue(ijk, newVal);
        }
        if (bboxAcc.isValueOn(ijm1k) || bboxAcc.isValueOn(ijp1k)) {
            auto val = vNextAcc.getValue(ijk);
            Vec3s newVal = Vec3s(val[0], 0, val[2]);
            vNextAcc.setValue(ijk, newVal);
        }
        if (bboxAcc.isValueOn(ijkm1) || bboxAcc.isValueOn(ijkp1)) {
            auto val = vNextAcc.getValue(ijk);
            Vec3s newVal = Vec3s(val[0], val[1], 0);
            vNextAcc.setValue(ijk, newVal);
        }
    }


    BoolTree::Ptr interiorMask(new BoolTree(false));
    interiorMask->topologyUnion(mVCurr->tree());
    tools::erodeActiveValues(*interiorMask, /*iterations=*/1, tools::NN_FACE, tools::IGNORE_TILES);
    BoolGrid::Ptr interiorGrid = BoolGrid::create(interiorMask);
    interiorGrid->setTransform(mXform);
    mInterior = interiorGrid->copy();
    mInterior->setName("interior");

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

    auto boolAcc = mInterior->getAccessor();
    int count = 0;
    for (auto iter = mVCurr->beginValueOn(); iter; ++iter) {
        math::Coord ijk = iter.getCoord();

        if (boolAcc.isValueOn(ijk)) {
            auto val = vNextAcc.getValue(ijk);
            std::cout << "vnext = " << ijk << " = " << val << std::endl;
        }
    }
    std::cout << "\t== divergence after vel bdry crct = " << divAfter << std::endl;
}


void
FlipSolver::pressureProjection() {
    using TreeType = FloatTree;
    using ValueType = TreeType::ValueType;

    const ValueType zero = zeroVal<ValueType>();
    const double epsilon = math::Delta<ValueType>::value();

    BoolTree::Ptr interiorMask(new BoolTree(false));
    interiorMask->topologyUnion(mVCurr->tree());
    tools::erodeActiveValues(*interiorMask, /*iterations=*/1, tools::NN_FACE, tools::IGNORE_TILES);
    BoolGrid::Ptr interiorGrid = BoolGrid::create(interiorMask);
    interiorGrid->setTransform(mXform);

    mDivBefore = tools::divergence(*mVCurr);
    mDivBefore->setName("div_before");
    (mDivBefore->tree()).topologyIntersection(interiorGrid->tree());
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

    math::pcg::State state = math::pcg::terminationDefaults<ValueType>();
    state.iterations = 100000;
    state.relativeError = state.absoluteError = epsilon;
    FlipSolver::BoundaryOp bop(mVoxelSize, mBBoxLS, mCollider, mVCurr);

    util::NullInterrupter interrupter;
    FloatTree::Ptr fluidPressure = tools::poisson::solveWithBoundaryConditions(
        mDivBefore->tree(), bop, state, interrupter, /*staggered=*/true);
    FloatGrid::Ptr fluidPressureGrid = FloatGrid::create(fluidPressure);
    fluidPressureGrid->setTransform(mXform);
    mPressure = fluidPressureGrid->copy();
    mPressure->setName("pressure");
    (fluidPressureGrid->tree()).topologyIntersection(interiorGrid->tree());

    Vec3SGrid::Ptr grad = tools::gradient(*fluidPressureGrid);
    grad->setGridClass(GRID_STAGGERED);
    (grad->tree()).topologyIntersection(interiorGrid->tree());

    auto vNextAcc = mVNext->getAccessor();
    auto vCurrAcc = mVCurr->getAccessor();
    auto gradAcc = grad->getAccessor();
    auto boolAcc = interiorGrid->getAccessor();

    for (auto iter = mVNext->beginValueOn(); iter; ++iter) {
        math::Coord ijk = iter.getCoord();

        auto val = vCurrAcc.getValue(ijk) - gradAcc.getValue(ijk) * mVoxelSize * mVoxelSize;
        //std::cout << "vCurr = " << vCurrAcc.getValue(ijk) << "\tgradAcc.getValue(ijk) = " << gradAcc.getValue(ijk) << std::endl;
        vNextAcc.setValue(ijk, val);
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
    std::cout << "\t== divergence after " << divAfter << std::endl;

    std::cout << "Success: " << state.success << std::endl;
    std::cout << "Iterations: " << state.iterations << std::endl;
    std::cout << "Relative error: " << state.relativeError << std::endl;
    std::cout << "Absolute error: " << state.absoluteError << std::endl;
    std::cout << "before dilate solution->activeVoxelCount() =  " << fluidPressure->activeVoxelCount() << std::endl;
}


void
FlipSolver::pressureProjection4() {
    std::cout << "pressure projection 4 begins" << std::endl;
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
    //std::cout << "\t== divergence before " << divBefore << std::endl;

    MaskGridType* domainMaskGrid = new MaskGridType(*mDivBefore); // match input grid's topology
    domainMaskGrid->topologyDifference(*mBBoxLS);
    // domainMaskGrid->topologyDifference(*mCollider);

    math::pcg::State state = math::pcg::terminationDefaults<ValueType>();
    state.iterations = 100000;
    state.relativeError = state.absoluteError = epsilon;
    FlipSolver::BoundaryOp bop(mVoxelSize, mBBoxLS, mCollider, mVCurr);

    util::NullInterrupter interrupter;
    // FloatTree::Ptr fluidPressure = tools::poisson::solveWithBoundaryConditions(
    //     mDivBefore->tree(), bop, state, interrupter, /*staggered=*/true);

    FloatTree::Ptr fluidPressure = tools::poisson::solveWithBoundaryConditionsAndPreconditioner<PCT>(
        mDivBefore->tree(), domainMaskGrid->tree(), bop, state, interrupter, /*staggered=*/true);
    //         divGrid->tree(), domainMaskGrid->tree(), boundaryOp, parms.outputState,
    //             *parms.interrupter, staggered);

    // FloatTree::Ptr fluidPressure = tools::poisson::solveWithBoundaryConditions(
    FloatGrid::Ptr fluidPressureGrid = FloatGrid::create(fluidPressure);
    fluidPressureGrid->setTransform(mXform);
    mPressure = fluidPressureGrid->copy();
    mPressure->setName("pressure");
    // (fluidPressureGrid->tree()).topologyIntersection(interiorGrid->tree());

    Vec3SGrid::Ptr grad = tools::gradient(*fluidPressureGrid);
    grad->setGridClass(GRID_STAGGERED);
    // (grad->tree()).topologyIntersection(interiorGrid->tree());
    // NOTE: line 712-714 in SOP_OpenVDB_Remove_Divergence
    grad->topologyUnion(*mVCurr);
    grad->topologyIntersection(*mVCurr);
    openvdb::tools::pruneInactive(grad->tree());

    auto vNextAcc = mVNext->getAccessor();
    auto vCurrAcc = mVCurr->getAccessor();
    auto gradAcc = grad->getAccessor();
    auto boolAcc = interiorGrid->getAccessor();

    int count = 0;
    for (auto iter = mVCurr->beginValueOn(); iter; ++iter) {
        math::Coord ijk = iter.getCoord();

        auto val = vCurrAcc.getValue(ijk) - gradAcc.getValue(ijk) * mVoxelSize * mVoxelSize;
        //std::cout << "vCurr = " << vCurrAcc.getValue(ijk) << "\tgradAcc.getValue(ijk) = " << gradAcc.getValue(ijk) << std::endl;
        if (boolAcc.isValueOn(ijk)) {
            vNextAcc.setValue(ijk, val);
            // std::cout << "newvel = " << ijk << " = " << val << "\tvcurr = " << vCurrAcc.getValue(ijk) << std::endl;
        }
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
    std::cout << "pressure projection 4 ends" << std::endl;
}


void
FlipSolver::pressureProjection5() {
    std::cout << "pressure projection 5 begins" << std::endl;
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
    // domainMaskGrid->topologyDifference(*mCollider);

    math::pcg::State state = math::pcg::terminationDefaults<ValueType>();
    state.iterations = 100000;
    state.relativeError = state.absoluteError = epsilon;
    FlipSolver::BoundaryOp bop(mVoxelSize, mBBoxLS, mCollider, mVCurr);

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
        std::cout << "pressure " << ijk << " = " << val << " div = " << divijk << std::endl;
    }

    // From conversation with Greg
    tools::dilateActiveValues(*fluidPressure, /*iterations=*/1, tools::NN_FACE, tools::IGNORE_TILES);
    mPressure = fluidPressureGrid->copy();
    mPressure->setName("pressure");
    // (fluidPressureGrid->tree()).topologyIntersection(interiorGrid->tree());

    GradientStaggered<FloatGrid> gradientOp(*fluidPressureGrid);
    Vec3SGrid::Ptr grad = gradientOp.process();
    grad->setGridClass(GRID_STAGGERED);
    
    auto gradAccOne = grad->getAccessor();
    auto vCurrAcc = mVCurr->getAccessor();
    // for (auto iter = grad->beginValueOn(); iter; ++iter) {
    //     math::Coord ijk = iter.getCoord();

    //     auto val = gradAccOne.getValue(ijk);
    //     std::cout << "grad acc one " << ijk << " = " << val * mVoxelSize << " vCurr = " << vCurrAcc.getValue(ijk) << std::endl;
    // }


    // (grad->tree()).topologyIntersection(interiorGrid->tree());
    // NOTE: line 712-714 in SOP_OpenVDB_Remove_Divergence
    grad->topologyUnion(*mVCurr);
    grad->topologyIntersection(*mVCurr);
    openvdb::tools::pruneInactive(grad->tree());

    auto vNextAcc = mVNext->getAccessor();
    auto gradAcc = grad->getAccessor();
    auto boolAcc = interiorGrid->getAccessor();

    int count = 0;
    for (auto iter = mVCurr->beginValueOn(); iter; ++iter) {
        math::Coord ijk = iter.getCoord();
        Vec3s gradijk;
        gradijk[0] = pressureAcc.getValue(ijk) - pressureAcc.getValue(ijk.offsetBy(-1, 0, 0));
        gradijk[1] = pressureAcc.getValue(ijk) - pressureAcc.getValue(ijk.offsetBy(0, -1, 0));
        gradijk[2] = pressureAcc.getValue(ijk) - pressureAcc.getValue(ijk.offsetBy(0, 0, -1));

        if (boolAcc.isValueOn(ijk)) {
            auto val = vCurrAcc.getValue(ijk) - gradijk * mVoxelSize;
            vNextAcc.setValue(ijk, val);
            // This is only multiplied by mVoxelSize because in the computation of
            // gradijk, I don't divide by mVoxelSize.
            std::cout << "gradijk = " << ijk << " = " << gradijk * mVoxelSize << "\tvnext = " << vNextAcc.getValue(ijk) << std::endl;
        }
        //std::cout << "vCurr = " << vCurrAcc.getValue(ijk) << "\tgradAcc.getValue(ijk) = " << gradAcc.getValue(ijk) << std::endl;
        // if (boolAcc.isValueOn(ijk)) {
        //     vNextAcc.setValue(ijk, val);
        //     std::cout << "newvel = " << ijk << " = " << val << "\tvcurr = " << vCurrAcc.getValue(ijk) << std::endl;
        // }
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
    std::cout << "pressure projection 5 ends" << std::endl;
}


void
FlipSolver::gridVelocityUpdate(float const dt) {
    addGravity(dt);
    velocityBCCorrection();
    pressureProjection5();
    velocityBCCorrection();
}


void
FlipSolver::substep(float const dt) {
    particlesToGrid();
    gridVelocityUpdate(dt);
    gridToParticles();
    advectParticles(dt);
}


void
FlipSolver::render() {
    float const dt = 1.f/24.f;
    for (int frame = 0; frame < 1; ++frame) {
        std::cout << "frame = " << frame << "\n";
        substep(dt);
        writeVDBs(frame);
        writeVDBsDebug(frame);
    }
}


void
FlipSolver::particlesToGrid2(){
    auto const xform =
        math::Transform::createLinearTransform(mVoxelSize);

    // Create a vector with four point positions.
    std::vector<Vec3s> positions;
    positions.push_back(Vec3s(0.f, 0.f, 0.f));
    positions.push_back(Vec3s(1.5f * 0.1, 0.f, 0.f));
    std::vector<Vec3s> velocities;
    velocities.push_back(Vec3s(1.f, 1.f, 1.f));
    velocities.push_back(Vec3s(1.f, 2.f, 3.f));

    points::PointAttributeVector<Vec3s> positionsWrapper(positions);


    // Create a PointDataGrid
    points::PointDataGrid::Ptr points =
        points::createPointDataGrid<points::NullCodec,
                        points::PointDataGrid>(positions, *xform);
    tools::PointIndexGrid::Ptr pointIndexGrid =
        openvdb::tools::createPointIndexGrid<tools::PointIndexGrid>(
            positionsWrapper, *xform);

    points::TypedAttributeArray<Vec3s, points::NullCodec>::registerType();
    NamePair velocityAttribute =
        openvdb::points::TypedAttributeArray<Vec3s, points::NullCodec>::attributeType();
    openvdb::points::appendAttribute(points->tree(), "velocity", velocityAttribute);
    points::PointAttributeVector<Vec3s> velocityWrapper(velocities);
    points::populateAttribute<points::PointDataTree,
        tools::PointIndexTree, points::PointAttributeVector<Vec3s>>(
            points->tree(), pointIndexGrid->tree(), "velocity", velocityWrapper);

    TreeBase::Ptr tree =
        points::rasterizeTrilinear</*staggered=*/true, Vec3s>(points->tree(), "velocity");

    Vec3STree::Ptr velTree = DynamicPtrCast<Vec3STree>(tree);
    mVCurr = Vec3SGrid::create(velTree);
    mVCurr->setGridClass(GRID_STAGGERED);
    mVCurr->setTransform(xform);
    auto vCurrAcc = mVCurr->getAccessor();

    for (auto iter = mVCurr->beginValueOn(); iter; ++iter) {
        auto ijk = iter.getCoord();
        std::cout << "=== mVCurr" << ijk << " = " << vCurrAcc.getValue(ijk) << std::endl;
    }
    mFlags = Int32Grid::create(0);
    (mFlags->tree()).topologyUnion(mVCurr->tree());
    mFlags->setTransform(xform);
    Int32Grid::Accessor flagAcc = mFlags->getAccessor();

    for (auto iter = mFlags->beginValueOn(); iter; ++iter) {
        math::Coord ijk = iter.getCoord();
        if (ijk.y() <= -1 || ijk.x() <= -1 || ijk.x() >= 3 || ijk.z() <= -1 || ijk.z() >= 1) {
            flagAcc.setValue(ijk, 2); // Neumann pressure
        }

        if (ijk.y() >= 1) {
            flagAcc.setValue(ijk, 1); // Dirichlet pressure
        }
    }
    // for (int k = -1; k <= 1; ++k) {
    //     for (int j = -1; j <= 1; ++j) {
    //         for (int i = -1; i <= 3; ++i) {
    //             math::Coord ijk(i, j, k);
    //             std::cout << "mFlags " << ijk << " = " << flagAcc.getValue(ijk) << std::endl; // Dirichlet pressure
    //         }
    //     }
    // }
}
    

void
FlipSolver::gridToParticles(){}


void
FlipSolver::pressureProjection2(){
    using TreeType = FloatTree;
    using ValueType = TreeType::ValueType;

    const ValueType zero = zeroVal<ValueType>();
    const double epsilon = math::Delta<ValueType>::value();


    std::cout << "\n\nflip::pressureProjection2 begins" << std::endl;
    auto const vCurrAcc = mVCurr->getConstAccessor();
    // for (int k = -1; k < 3; ++k) {
    //     for (int j = -1; j < 3; ++j) {
    //         for (int i = -1; i < 3; ++i) {
    //             math::Coord ijk(i, j, k);
    //             auto val = vCurrAcc.getValue(ijk);
    //             if (val.length() > 1.0e-5) {
    //                 std::cout << "vel" << ijk << " = " << vCurrAcc.getValue(ijk) << std::endl;
    //             }
    //         }
    //     }
    // }

    BoolTree::Ptr interiorMask(new BoolTree(false));
    interiorMask->topologyUnion(mVCurr->tree());
    // This is one way to do this, another way is to deduce the topology from voxels
    // where the velocity is not-zero, which might be more accurate.
    tools::erodeActiveValues(*interiorMask, /*iterations=*/1, tools::NN_FACE, tools::IGNORE_TILES);
    BoolGrid::Ptr interiorGrid = BoolGrid::create(interiorMask);
    interiorGrid->setTransform(mVCurr->transform().copy());
    BoolGrid::ConstAccessor intrAcc = interiorGrid->getConstAccessor();

    // for (auto iter = interiorGrid->beginValueOn(); iter; ++iter) {
    //     math::Coord ijk = iter.getCoord();
    //     auto val = intrAcc.getValue(ijk);
    //     std::cout << "interior " << ijk << " = " << val << std::endl;
    // }

    // Setup the right hand side 
    FloatGrid::Ptr divGrid = tools::divergence(*mVCurr);
    //(divGrid->tree()).topologyIntersection(interiorGrid->tree());
    auto divAcc = divGrid->getConstAccessor();
    std::cout << "divgrid before pressure projection" << std::endl;
    // Note that ijk = [0,0,0] is not printed because the divergence in that cell is 0
    for (auto iter = divGrid->beginValueOn(); iter; ++iter) {
        math::Coord ijk = iter.getCoord();
        auto val = divAcc.getValue(ijk);
        std::cout << "beginvalon div " << ijk << " = " << val << std::endl;
    }

    // Setup conjugate gradient
    std::cout << "RHS: divGrid->activeVoxelCount() = " << divGrid->activeVoxelCount() << std::endl;
    math::pcg::State state = math::pcg::terminationDefaults<ValueType>();
    state.iterations = 100000;
    state.relativeError = state.absoluteError = epsilon;
    util::NullInterrupter interrupter;
    FloatTree::Ptr fluidPressure = tools::poisson::solveWithBoundaryConditions(
        divGrid->tree(), FlipSolver::BoundaryFooBarOp(mVoxelSize), state, interrupter, /*staggered=*/true);
    FloatGrid::Ptr fluidPressureGrid = FloatGrid::create(fluidPressure);

    auto fluidPressureAcc = fluidPressureGrid->getAccessor();
    for (auto iter = fluidPressure->beginValueOn(); iter; ++iter) {
        math::Coord ijk = iter.getCoord();
        auto val = fluidPressureAcc.getValue(ijk);
        std::cout << "fluidPressure " << ijk << " = " << val << std::endl;
    }

    std::cout << "Success: " << state.success << std::endl;
    std::cout << "Iterations: " << state.iterations << std::endl;
    std::cout << "Relative error: " << state.relativeError << std::endl;
    std::cout << "Absolute error: " << state.absoluteError << std::endl;
    std::cout << "before dilate solution->activeVoxelCount() =  " << fluidPressure->activeVoxelCount() << std::endl;

    tools::dilateActiveValues(*fluidPressure, /*iterations=*/1, tools::NN_FACE, tools::IGNORE_TILES);
    
    fluidPressureGrid->setTransform(mVCurr->transform().copy());
    auto fluidAcc = fluidPressureGrid->getAccessor();
    math::Coord ijk(0, 0, 0);
    float p000 = 0.f;//fluidPressureAcc.getValue(ijk);
    ijk = math::Coord(1, 0, 0);
    float p100 = 0.f; //fluidPressureAcc.getValue(ijk);
    ijk = math::Coord(2, 0, 0);
    float p200 = 0.f; //fluidPressureAcc.getValue(ijk);

    ijk = math::Coord(-1, 0, 0);
    fluidAcc.setValue(ijk, p000);
    ijk = math::Coord(0, -1, 0);
    fluidAcc.setValue(ijk, p000);
    ijk = math::Coord(1, -1, 0);
    fluidAcc.setValue(ijk, p100);
    ijk = math::Coord(2, -1, 0);
    fluidAcc.setValue(ijk, p200);
    ijk = math::Coord(0, 0, -1);
    fluidAcc.setValue(ijk, p000);
    ijk = math::Coord(1, 0, -1);
    fluidAcc.setValue(ijk, p100);
    ijk = math::Coord(2, 0, -1);
    fluidAcc.setValue(ijk, p200);
    ijk = math::Coord(0, 0, 1);
    fluidAcc.setValue(ijk, p000);
    ijk = math::Coord(0, 1, 0);
    fluidAcc.setValue(ijk, 0);
    ijk = math::Coord(1, -1, 0);
    fluidAcc.setValue(ijk, p100);
    ijk = math::Coord(1, 0, 1);
    fluidAcc.setValue(ijk, p100);
    ijk = math::Coord(1, 1, 0);
    fluidAcc.setValue(ijk, 0);
    ijk = math::Coord(2, 0, 1);
    fluidAcc.setValue(ijk, p200);
    ijk = math::Coord(2, 1, 0);
    fluidAcc.setValue(ijk, 0);
    ijk = math::Coord(3, 0, 0);
    fluidAcc.setValue(ijk, p200);
    for (auto iter = fluidPressure->beginValueOn(); iter; ++iter) {
        math::Coord ijk = iter.getCoord();
        std::cout << "p" << ijk << " = " << fluidAcc.getValue(ijk) << std::endl;
    }
    // p[0, 0, 0] = 1.25
    // p[1, 0, 0] = 2.5
    // p[2, 0, 0] = 6.25

    Vec3SGrid::Ptr grad = tools::gradient(*fluidPressureGrid);
    grad->setGridClass(GRID_STAGGERED);
    // std::cout << "grad->activeVoxelCount() = " << grad->activeVoxelCount() << std::endl;
    std::cout << "grad->voxelSize() = " << grad->voxelSize() << std::endl;

    mVNext = mVCurr->copy();
    mVNext->setGridClass(GRID_STAGGERED);
    mVNext->setTransform(mVCurr->transform().copy());
    auto vNextAcc = mVNext->getAccessor();
    auto gradAcc = grad->getAccessor();

    for (auto iter = mVNext->beginValueOn(); iter; ++iter) {
        math::Coord ijk = iter.getCoord();
        auto val = vCurrAcc.getValue(ijk) - mVoxelSize * mVoxelSize * gradAcc.getValue(ijk);
        vNextAcc.setValue(ijk, val);
    }

    FloatGrid::Ptr divGridAfter = tools::divergence(*mVNext);
    (divGridAfter->tree()).topologyIntersection(interiorGrid->tree());
    FloatGrid::ConstAccessor divAfterAcc = divGridAfter->getConstAccessor();
    std::cout << "divgrid after pressure projection" << std::endl;
    // Note that ijk = [0,0,0] is not printed because the divergence in that cell is 0
    for (auto iter = divGridAfter->beginValueOn(); iter; ++iter) {
        math::Coord ijk = iter.getCoord();
        auto val = divAfterAcc.getValue(ijk);
        std::cout << "beginvalon div after " << ijk << " = " << val << std::endl;
    }

    // //for (auto iter = mVCurr->beginValueOn(); iter; ++iter) {
    // //math::Coord ijk = iter.getCoord();
    // std::cout << "=== Debug mVoxelSize = " << mVoxelSize <<
    //     " mVNext->voxelSize() = " << mVNext->voxelSize() <<
    //     " grad->voxelSize() = " << grad->voxelSize() <<
    //     std::endl;
    // for (int k = 0; k <=1; ++k) {
    // for (int j = 0; j <= 1; ++j) {
    // for (int i = 0; i <= 2; ++i) {
    //     math::Coord ijk(i, j, k);
    //     auto val = vNextAcc.getValue(ijk) - gradAcc.getValue(ijk) * mVoxelSize /** mVoxelSize */;
    //     auto gradVal = gradAcc.getValue(ijk);
    //     auto scaleGrad1 = gradVal * mVoxelSize;
    //     auto scaleGrad2 = gradVal * mVoxelSize * mVoxelSize;
    //     std::cout << "vNext.getValue(ijk) = " << vNextAcc.getValue(ijk) << std::endl;
    //     std::cout << "gradVal = " << gradVal << 
    //         " scaleGrad1 = " << scaleGrad1 <<
    //         " scaleGrad2 = " << scaleGrad2 <<
    //         std::endl;
    //     vNextAcc.setValue(ijk, val);
    // }
    // }
    // }

    // // Div grid after
    // FloatGrid::Ptr divGridAfter = tools::divergence(*mVNext);
    // (divGridAfter->tree()).topologyIntersection(interiorGrid->tree());
    // FloatGrid::ConstAccessor divAfterAcc = divGridAfter->getConstAccessor();
    // std::cout << "divgrid after pressure projection" << std::endl;
    // // Note that ijk = [0,0,0] is not printed because the divergence in that cell is 0
    // for (auto iter = divGridAfter->beginValueOn(); iter; ++iter) {
    //     math::Coord ijk = iter.getCoord();
    //     auto val = divAfterAcc.getValue(ijk);
    //     std::cout << "beginvalon div after " << ijk << " = " << val << std::endl;
    // }


    // // Make a deep copy of fluidPressure
    // TreeBase::Ptr baseFullPressure = fluidPressure->copy();
    // // Cast down to the concrete type to query values
    // openvdb::FloatTree *fullPressure = dynamic_cast<FloatTree*>(baseFullPressure.get());
    // tools::dilateActiveValues(*fullPressure, /*iterations=*/1, tools::NN_FACE, tools::IGNORE_TILES);

    // // TODO: parallelize this:

    // std::cout << "flip::pressureProjection2 ends" << std::endl;
}

void
FlipSolver::pressureProjection3(){
    std::cout << "flip::pressureProjection begin" << std::endl;
    BoolTree::Ptr interiorMask(new BoolTree(false));
    interiorMask->topologyUnion(mVCurr->tree());
    tools::erodeActiveValues(*interiorMask, /*iterations=*/1, tools::NN_FACE, tools::IGNORE_TILES);
    BoolGrid::Ptr interiorGrid = BoolGrid::create(interiorMask);
    BoolGrid::ConstAccessor intrAcc = interiorGrid->getConstAccessor();

    for (auto iter = interiorGrid->beginValueOn(); iter; ++iter) {
        math::Coord ijk = iter.getCoord();
        auto val = intrAcc.getValue(ijk);
        std::cout << "ijk = " << ijk << ", val = " << val << std::endl;
    }

    FloatGrid::Ptr divGrid = tools::divergence(*mVCurr);
    (divGrid->tree()).topologyIntersection(interiorGrid->tree());
    FloatGrid::ConstAccessor divAcc = divGrid->getConstAccessor();
    std::cout << "divgrid before pressure projection" << std::endl;
    // Note that ijk = [0,0, 0] is not printed because the divergence in that cell is 0
    for (auto iter = divGrid->beginValueOn(); iter; ++iter) {
        math::Coord ijk = iter.getCoord();
        auto val = divAcc.getValue(ijk);
        std::cout << "div ijk = " << ijk << ", val = " << val << std::endl;
    }

    BoolTree::Ptr pressureMask(new BoolTree(false));
    pressureMask->topologyUnion(divGrid->tree());
    tools::dilateActiveValues(*pressureMask, /*iterations=*/1, tools::NN_FACE, tools::IGNORE_TILES);
    BoolGrid::Ptr pressureGridMask = BoolGrid::create(pressureMask);
    BoolGrid::ConstAccessor prsAcc = pressureGridMask->getConstAccessor();
    std::cout << "pressuregrid after pressure projection" << std::endl;
    for (auto iter = pressureGridMask->beginValueOn(); iter; ++iter) {
        math::Coord ijk = iter.getCoord();
        auto val = prsAcc.getValue(ijk);
        std::cout << "pres ijk = " << ijk << ", val = " << val << std::endl;
    }



    std::cout << "flip::pressureProjection end" << std::endl;

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

void openvdb_points_for_houndstooth() {
    std::vector<std::string> vdb_names = {"waterfall_100k.vdb",
                                          "waterfall_1mil.vdb",
                                          "waterfall_10mil.vdb"};

    std::vector<std::string> output_names = {"waterfall_100k.txt",
                                             "waterfall_1mil.txt",
                                             "waterfall_10mil.txt"};

    std::vector<int> point_numbers = {100000, 1000000, 10000000};

    int const LENGTH = 3;

    // WRITE TO FILES
    for (int i = 0; i < LENGTH; ++i) {

        auto const vdb_name = vdb_names[i];
        auto const output_name = output_names[i];

        std::cout << "Writing " << vdb_name << std::endl;

        std::ofstream outputFile(output_name.c_str());

        openvdb::io::File file(vdb_name.c_str());
        // Open the file. This reads the file header, but not any grids.
        file.open();
        // Read the grid by name.
        openvdb::GridBase::Ptr baseGrid = file.readGrid("points");
        file.close();
        // From the example above, "points" is known to be a PointDataGrid,
        // so cast the generic grid pointer to a PointDataGrid pointer.
        auto grid = openvdb::gridPtrCast<openvdb::points::PointDataGrid>(baseGrid);
        openvdb::Index64 count = openvdb::points::pointCount(grid->tree());
        std::cout << "PointCount=" << count << std::endl;
        // Iterate over all the leaf nodes in the grid.
        for (auto leafIter = grid->tree().cbeginLeaf(); leafIter; ++leafIter) {
            // Verify the leaf origin.
            // Extract the position attribute from the leaf by name (P is position).
            const openvdb::points::AttributeArray& array =
                leafIter->constAttributeArray("P");
            // Create a read-only AttributeHandle. Position always uses Vec3f.
            openvdb::points::AttributeHandle<openvdb::Vec3f> positionHandle(array);
            // Iterate over the point indices in the leaf.
            for (auto indexIter = leafIter->beginIndexOn(); indexIter; ++indexIter) {
                // Extract the voxel-space position of the point.
                openvdb::Vec3f voxelPosition = positionHandle.get(*indexIter);
                // Extract the index-space position of the voxel.
                const openvdb::Vec3d xyz = indexIter.getCoord().asVec3d();
                // Compute the world-space position of the point.
                openvdb::Vec3f worldPosition =
                    grid->transform().indexToWorld(voxelPosition + xyz);
                // Verify the index and world-space position of the point
                outputFile << worldPosition[0] << "," << worldPosition[1] << "," << worldPosition[2] << "\n";
            }
        } // leafIter
        outputFile.close();
    } // LENGTH

    for (int i = 0; i < LENGTH; ++i) {
        std::string file_name = output_names[i];
        std::ifstream inputFile(file_name.c_str());
        std::string line;
        std::vector<Vector3> vectorList;

        if (inputFile.is_open()) {
            while (std::getline(inputFile, line)) {
                std::stringstream ss(line);
                std::string value;
                std::vector<float> values;

                while (std::getline(ss, value, ',')) {
                    values.push_back(std::stof(value));
                }

                if (values.size() == 3) {
                    Vector3 vector(values[0], values[1], values[2]);
                    vectorList.push_back(vector);
                }
            }

            inputFile.close();
        }
        std::cout << "i = " << i << " vector length = " << vectorList.size() << std::endl;
    }
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
