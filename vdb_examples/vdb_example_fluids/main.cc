// Copyright Contributors to the OpenVDB Project
// SPDX-License-Identifier: MPL-2.0

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <math.h>
#include <string>
#include <vector>

#include <openvdb/openvdb.h>
#include <openvdb/points/PointAdvect.h> // for advectPoints
#include <openvdb/points/PointAttribute.h> // for appendAttribute
#include <openvdb/points/PointConversion.h>
#include <openvdb/points/PointDataGrid.h> // for PointDataGrid
#include <openvdb/points/PointRasterizeTrilinear.h> // for rasterizing to the grid
#include <openvdb/points/PointSample.h> // for PointSample
#include <openvdb/points/PointScatter.h> // for point sampling
#include <openvdb/tools/Composite.h> // for tools::compMax
#include <openvdb/tools/GridOperators.h> // for divergence and gradient
#include <openvdb/tools/Interpolation.h> // for BoxSampler
#include <openvdb/tools/MeshToVolume.h> // for createLevelSetBox
#include <openvdb/tools/Morphology.h> // for erodeActiveValues
#include <openvdb/tools/PoissonSolver.h> // for poisson solve
#include <openvdb/tools/VolumeAdvect.h> // for tools::VolumeAdvection
#include <openvdb/tree/LeafManager.h> // for LeafManager
#include <openvdb/tree/NodeManager.h> // for post processing bool grid
#include <openvdb/util/logging.h>

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
    // back to the particle.
    void gridToParticles();
    void updateParticles(float const dt);
    void updateParticlesVelocity();

    // Update particle position based on velocity on the grid
    void advectParticles(float const dt);
    void projectParticlesOutOfCollider();
    void enforceParticleColliderVelocity();

    // Make the velocity on the grid to be divergence free
    void pressureProjection(bool print);

    void gridVelocityUpdate(float const dt);

    void velocityBCCorrection(Vec3SGrid& vecGrid);
    void extrapolateVelocity(Vec3SGrid& vecGrid, int const iterations);

    void addGravity(float const dt);
    void computeFlipVelocity();

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
            if (isInsideCollider) {
                double delta = 0.0;
                // Neumann pressure from bbox
                if (neighbor.x() + 1 == ijk.x() /* left x-face */) {
                    delta += /* voxelSize * */ vCurr->tree().getValue(ijk)[0];
                }
                if (neighbor.x() - 1 == ijk.x() /* right x-face */) {
                    delta -= /* voxelSize * */ vCurr->tree().getValue(neighbor)[0];
                }
                if (neighbor.y() + 1 == ijk.y() /* bottom y-face */) {
                    delta += /* voxelSize * */ vCurr->tree().getValue(ijk)[1];
                }
                if (neighbor.y() - 1 == ijk.y() /* top y-face */) {
                    delta -= /* voxelSize * */ vCurr->tree().getValue(neighbor)[1];
                }
                if (neighbor.z() + 1 == ijk.z() /* back z-face */) {
                    delta += /* voxelSize *  */ vCurr->tree().getValue(ijk)[2];
                }
                if (neighbor.z() - 1 == ijk.z() /* front z-face */) {
                    delta -= /* voxelSize *  */ vCurr->tree().getValue(neighbor)[2];
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
                    const auto oldVel = velHandle.get(*indexIter);
                    const auto vPic = vPicHandle.get(*indexIter);
                    const auto vFlip = oldVel + vFlipHandle.get(*indexIter);
                    const auto newVel = alpha * vPic + (1 - alpha) * vFlip;
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
    // vOld. To be used with foreach in LeafManager.
    struct ComputeFlipVelocityOp
    {
        ComputeFlipVelocityOp(Vec3SGrid::Ptr vOld,
                              Vec3SGrid::Ptr vNext) :
                              vOld(vOld),
                              vNext(vNext) {}

        template <typename T>
        void operator()(T &leaf, size_t) const
        {
            auto vOldAcc = vOld->getAccessor();
            auto vNextAcc = vNext->getAccessor();
            for (typename T::ValueOnIter iter = leaf.beginValueOn(); iter; ++iter) {
                auto ijk = iter.getCoord();
                Vec3s val = vNextAcc.getValue(ijk) - vOldAcc.getValue(ijk);
                iter.setValue(val);
            }
        }

        Vec3SGrid::Ptr vOld;
        Vec3SGrid::Ptr vNext;
    };// ComputeFlipVelocityOp


    struct ColliderSDFSampler
    {
        ColliderSDFSampler(FloatGrid::ConstPtr collider) : collider(collider) {}

        double sample(Vec3d const& position) const
        {
            return tools::BoxSampler::sample(
                collider->tree(), collider->transform().worldToIndex(position));
        }

        Vec3d normal(Vec3d const& position) const
        {
            const Vec3d ijk = collider->transform().worldToIndex(position);
            Vec3d grad(
                tools::BoxSampler::sample(collider->tree(), ijk + Vec3d(1, 0, 0)) -
                    tools::BoxSampler::sample(collider->tree(), ijk - Vec3d(1, 0, 0)),
                tools::BoxSampler::sample(collider->tree(), ijk + Vec3d(0, 1, 0)) -
                    tools::BoxSampler::sample(collider->tree(), ijk - Vec3d(0, 1, 0)),
                tools::BoxSampler::sample(collider->tree(), ijk + Vec3d(0, 0, 1)) -
                    tools::BoxSampler::sample(collider->tree(), ijk - Vec3d(0, 0, 1)));
            const double length = grad.length();
            if (length > 1.0e-12) return grad / length;
            return Vec3d(0.0, 1.0, 0.0);
        }

        FloatGrid::ConstPtr collider;
    };


    struct ProjectParticleOutOfColliderOp
    {
        ProjectParticleOutOfColliderOp(FloatGrid::ConstPtr collider,
                                       double const padding) :
                                       sampler(collider),
                                       padding(padding) {}

        template <typename T>
        void reset(T&, size_t) {}

        template <typename IndexIterT>
        void apply(Vec3d& position, const IndexIterT&) const
        {
            const double phi = sampler.sample(position);
            if (phi < padding) {
                position += (padding - phi) * sampler.normal(position);
            }
        }

        ColliderSDFSampler sampler;
        double padding;
    };


    struct EnforceParticleColliderVelocityOp
    {
        EnforceParticleColliderVelocityOp(Index64 const posAtrIdx,
                                          Index64 const velAtrIdx,
                                          math::Transform const& transform,
                                          FloatGrid::ConstPtr collider,
                                          double const padding) :
                                          posAtrIdx(posAtrIdx),
                                          velAtrIdx(velAtrIdx),
                                          transform(transform),
                                          sampler(collider),
                                          padding(padding) {}

        void operator()(const tree::LeafManager<points::PointDataTree>::LeafRange& range) const {
            for (auto leafIter = range.begin(); leafIter; ++leafIter) {
                points::AttributeArray const& posArray = leafIter->constAttributeArray(posAtrIdx);
                points::AttributeArray& velArray = leafIter->attributeArray(velAtrIdx);
                points::AttributeHandle<Vec3f> posHandle(posArray);
                points::AttributeWriteHandle<Vec3s> velHandle(velArray);

                for (auto indexIter = leafIter->beginIndexOn(); indexIter; ++indexIter) {
                    const Vec3d indexPosition =
                        indexIter.getCoord().asVec3d() + Vec3d(posHandle.get(*indexIter));
                    const Vec3d worldPosition = transform.indexToWorld(indexPosition);
                    Vec3s velocity = velHandle.get(*indexIter);
                    const double phi = sampler.sample(worldPosition);

                    if (phi < padding) {
                        const Vec3d n = sampler.normal(worldPosition);
                        const double normalVelocity = Vec3d(velocity).dot(n);
                        if (normalVelocity < 0.0) {
                            velocity = Vec3s(Vec3d(velocity) - normalVelocity * n);
                        }
                    }
                    velHandle.set(*indexIter, velocity);
                }
            }
        }

        Index64 posAtrIdx;
        Index64 velAtrIdx;
        math::Transform const& transform;
        ColliderSDFSampler sampler;
        double padding;
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
    Vec3SGrid::Ptr mVOld;
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

    openvdb::Index64 count = openvdb::points::pointCount(mPoints->tree());
    std::cout << "PointCount=" << count << std::endl;
}


void
FlipSolver::initializeDamBreak() {
    using BBox = math::BBox<Vec3s>;

    mXform = math::Transform::createLinearTransform(mVoxelSize);
    float const padding = 2 * mVoxelSize;

    Vec3s minFI = Vec3s(mVoxelSize, mVoxelSize, mVoxelSize);
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

    Vec3s obstacleMin = Vec3s(1.1f, 0.1f, 1.0f);
    Vec3s obstacleMax = Vec3s(3.1f, 4.1f, 4.1f);
    Coord obstacleMinCoord = mXform->worldToIndexNodeCentered(obstacleMin);
    Coord obstacleMaxCoord = mXform->worldToIndexNodeCentered(obstacleMax);
    FloatGrid::Ptr obstacleMask = FloatGrid::create(/*bg = */0.f);
    obstacleMask->denseFill(
        CoordBBox(obstacleMinCoord, obstacleMaxCoord), /*value = */ 1.0, /*active = */ true);
    obstacleMask->setTransform(mXform);

    Vec3s minBBoxvec = Vec3s(-padding, -padding, -padding);
    Vec3s maxBBoxvec = Vec3s(14.f + padding + mVoxelSize, 5.f + padding + mVoxelSize, 5.f + padding + mVoxelSize);
    Coord minBBoxcoord = mXform->worldToIndexNodeCentered(minBBoxvec);
    Coord maxBBoxcoord = mXform->worldToIndexNodeCentered(maxBBoxvec);
    mBBoxLS = FloatGrid::create(/*bg = */0.f);
    mBBoxLS->denseFill(CoordBBox(minBBoxcoord, maxBBoxcoord), /*value = */ 1.0, /*active = */ true);
    mBBoxLS->setTransform(mXform);
    mBBoxLS->topologyDifference(*negativeSpace);
    mBBoxLS->topologyDifference(*fluidLSInit);
    mBBoxLS->topologyUnion(*obstacleMask);
    mBBoxLS->setName("collider");
    openvdb::tools::pruneInactive(mBBoxLS->tree());

    const Vec3d domainMin = mXform->indexToWorld(minFICoord.asVec3d() - Vec3d(0.5));
    const Vec3d domainMax = mXform->indexToWorld(maxFIIntrCoord.asVec3d() + Vec3d(0.5));
    const Vec3d obstacleDomainMin = mXform->indexToWorld(obstacleMinCoord.asVec3d() - Vec3d(0.5));
    const Vec3d obstacleDomainMax = mXform->indexToWorld(obstacleMaxCoord.asVec3d() + Vec3d(0.5));
    mCollider = FloatGrid::create(-padding);
    mCollider->setTransform(mXform);
    mCollider->setGridClass(GRID_LEVEL_SET);
    mCollider->setName("collider_sdf");
    auto colliderAcc = mCollider->getAccessor();
    auto boxInteriorPhi = [](Vec3d const& p, Vec3d const& boxMin, Vec3d const& boxMax) {
        double outsideDistance2 = 0.0;
        double insideDistance = std::numeric_limits<double>::max();
        for (int axis = 0; axis < 3; ++axis) {
            if (p[axis] < boxMin[axis]) {
                const double d = boxMin[axis] - p[axis];
                outsideDistance2 += d * d;
            } else if (p[axis] > boxMax[axis]) {
                const double d = p[axis] - boxMax[axis];
                outsideDistance2 += d * d;
            } else {
                insideDistance = std::min(
                    insideDistance,
                    std::min(p[axis] - boxMin[axis], boxMax[axis] - p[axis]));
            }
        }
        return outsideDistance2 > 0.0 ? -std::sqrt(outsideDistance2) : insideDistance;
    };
    for (CoordBBox::Iterator<true> iter(CoordBBox(minBBoxcoord, maxBBoxcoord)); iter; ++iter) {
        const Vec3d p = mXform->indexToWorld((*iter).asVec3d());
        const double containerPhi = boxInteriorPhi(p, domainMin, domainMax);
        const double obstaclePhi = -boxInteriorPhi(p, obstacleDomainMin, obstacleDomainMax);
        const double phi = std::min(containerPhi, obstaclePhi);
        colliderAcc.setValue(*iter, static_cast<float>(phi));
    }

    FloatGrid::Ptr fluidScatterMask = fluidLSInit->deepCopy();
    fluidScatterMask->topologyDifference(*mBBoxLS);
    openvdb::tools::pruneInactive(fluidScatterMask->tree());

    mPoints = points::denseUniformPointScatter(*fluidScatterMask, mPointsPerVoxel);
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

    mVOld = mVCurr->deepCopy();
    mVOld->setName("v_old");

    mVNext = Vec3SGrid::create(Vec3s(0.f, 0.f, 0.f));
    (mVNext->tree()).topologyUnion(mVCurr->tree());
    mVNext->setGridClass(GRID_STAGGERED);
    mVNext->setTransform(mXform);
    mVNext->setName("v_next");
}


void
FlipSolver::addGravity(float const dt) {
    tree::LeafManager<Vec3STree> r(mVCurr->tree());
    FlipSolver::ApplyGravityOp op(dt, mGravity);
    r.foreach(op);
}


void
FlipSolver::computeFlipVelocity() {
    mVDiff = Vec3SGrid::create(Vec3s(0.f, 0.f, 0.f));
    (mVDiff->tree()).topologyUnion(mVCurr->tree());
    mVDiff->setGridClass(GRID_STAGGERED);
    mVDiff->setTransform(mXform);

    tree::LeafManager<Vec3STree> r(mVDiff->tree());
    FlipSolver::ComputeFlipVelocityOp op(mVOld, mVNext);
    r.foreach(op);
}


void
FlipSolver::velocityBCCorrection(Vec3SGrid& vecGrid) {
    auto acc = vecGrid.getAccessor();
    auto bboxAcc = mBBoxLS->getAccessor();

    for (auto iter = vecGrid.beginValueOn(); iter; ++iter) {
        math::Coord ijk = iter.getCoord();
        math::Coord im1jk = ijk.offsetBy(-1, 0, 0);
        math::Coord ijm1k = ijk.offsetBy(0, -1, 0);
        math::Coord ijkm1 = ijk.offsetBy(0, 0, -1);

        if (bboxAcc.isValueOn(im1jk) || bboxAcc.isValueOn(ijk)) {
            auto val = acc.getValue(ijk);
            Vec3s newVal = Vec3s(0, val[1], val[2]);
            acc.setValue(ijk, newVal);
        }
        if (bboxAcc.isValueOn(ijm1k) || bboxAcc.isValueOn(ijk)) {
            auto val = acc.getValue(ijk);
            Vec3s newVal = Vec3s(val[0], 0, val[2]);
            acc.setValue(ijk, newVal);
        }
        if (bboxAcc.isValueOn(ijkm1) || bboxAcc.isValueOn(ijk)) {
            auto val = acc.getValue(ijk);
            Vec3s newVal = Vec3s(val[0], val[1], 0);
            acc.setValue(ijk, newVal);
        }
    }
}


void
FlipSolver::extrapolateVelocity(Vec3SGrid& vecGrid, int const iterations) {
    if (iterations <= 0) return;

    const Coord offsets[6] = {
        Coord(-1, 0, 0), Coord(1, 0, 0),
        Coord(0, -1, 0), Coord(0, 1, 0),
        Coord(0, 0, -1), Coord(0, 0, 1)
    };

    for (int iteration = 0; iteration < iterations; ++iteration) {
        Vec3SGrid::Ptr oldGrid = vecGrid.deepCopy();
        auto oldAcc = oldGrid->getConstAccessor();
        auto acc = vecGrid.getAccessor();
        ColliderSDFSampler colliderSampler(mCollider);
        auto isDeepInsideCollider = [this, &colliderSampler](Coord const& coord) {
            if (!mCollider) return false;
            const Vec3d position = mXform->indexToWorld(coord.asVec3d());
            return colliderSampler.sample(position) < -0.5 * double(mVoxelSize);
        };

        for (auto iter = oldGrid->beginValueOn(); iter; ++iter) {
            const Coord ijk = iter.getCoord();

            for (Coord const& offset : offsets) {
                const Coord neighbor = ijk + offset;
                if (oldAcc.isValueOn(neighbor)) continue;
                if (isDeepInsideCollider(neighbor)) continue;

                Vec3s sum(0.0f);
                int count = 0;
                for (Coord const& averageOffset : offsets) {
                    const Coord averageCoord = neighbor + averageOffset;
                    if (isDeepInsideCollider(averageCoord)) continue;
                    if (oldAcc.isValueOn(averageCoord)) {
                        sum += oldAcc.getValue(averageCoord);
                        ++count;
                    }
                }
                if (count > 0) {
                    acc.setValue(neighbor, sum / float(count));
                }
            }
        }
    }
}


void
FlipSolver::pressureProjection(bool print) {
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
    mInterior = interiorGrid->copy();
    mInterior->setName("interior");

    mDivBefore = tools::divergence(*mVCurr);
    mDivBefore->setName("div_before");

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
    // Note: need to dilate in order to do one-sided difference
    // because we use a staggered grid velocity field.
    tools::dilateActiveValues(*fluidPressure, /*iterations=*/1, tools::NN_FACE, tools::IGNORE_TILES);

    fluidPressureGrid->setTransform(mXform);
    mPressure = fluidPressureGrid->copy();
    mPressure->setName("pressure");

    auto vCurrAcc = mVCurr->getAccessor();
    auto vNextAcc = mVNext->getAccessor();
    auto boolAcc = interiorGrid->getAccessor();
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
FlipSolver::gridVelocityUpdate(float const dt) {
    addGravity(dt);
    velocityBCCorrection(*mVCurr);
    pressureProjection(false /* print */);
    velocityBCCorrection(*mVNext);
    computeFlipVelocity();
    extrapolateVelocity(*mVNext, 3);
    velocityBCCorrection(*mVNext);
    extrapolateVelocity(*mVDiff, 3);
    velocityBCCorrection(*mVDiff);
    mDivAfter = tools::divergence(*mVNext);
    mDivAfter->setName("div_after");
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
    FlipSolver::FlipUpdateOp op(velIdx, vPicIdx, vFlipIdx, 0.05 /* PIC blend */);
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
        std::cout << "\nframe = " << frame << "\n";
        substep(dt);
        writeVDBs(frame);
        writeVDBsVerbose(frame);
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
    projectParticlesOutOfCollider();
}


void
FlipSolver::projectParticlesOutOfCollider() {
    FlipSolver::ProjectParticleOutOfColliderOp op(mCollider, 1.0e-4 * mVoxelSize);
    points::movePoints(*mPoints, op);
    enforceParticleColliderVelocity();
}


void
FlipSolver::enforceParticleColliderVelocity() {
    auto leafIter = (mPoints->tree()).beginLeaf();
    if (!leafIter) return;

    auto descriptor = leafIter->attributeSet().descriptor();
    Index64 posIdx = descriptor.find("P");
    Index64 velIdx = descriptor.find("velocity");

    tree::LeafManager<points::PointDataTree> leafManager(mPoints->tree());
    FlipSolver::EnforceParticleColliderVelocityOp op(
        posIdx, velIdx, *mXform, mCollider, 1.0e-3 * mVoxelSize);
    tbb::parallel_for(leafManager.leafRange(), op);
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
FlipSolver::writeVDBsVerbose(int const frame) {
    std::ostringstream ss;
    ss << "water_volume_" << std::setw(3) << std::setfill('0') << frame << ".vdb";
    std::string fileName(ss.str());
    openvdb::io::File file(fileName.c_str());

    openvdb::GridPtrVec grids;
    auto addGrid = [&grids](GridBase::Ptr const& grid) {
        if (grid) grids.push_back(grid);
    };
    addGrid(mPoints);
    addGrid(mBBoxLS);
    addGrid(mCollider);
    addGrid(mVOld);
    addGrid(mVCurr);
    addGrid(mVNext);
    addGrid(mVDiff);
    addGrid(mDivBefore);
    addGrid(mDivAfter);
    addGrid(mPressure);
    addGrid(mInterior);
    file.write(grids);
    file.close();
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

    FlipSolver flipSim(0.1f /* voxel size */);
    flipSim.render();
}
