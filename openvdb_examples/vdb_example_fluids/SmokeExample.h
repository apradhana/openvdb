// Copyright Contributors to the OpenVDB Project
// SPDX-License-Identifier: MPL-2.0
// 0 Neumann
// 1 interior
// 4 dirichlet
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
    void initialize2();
    void initialize4();

    void substep(float const dt, int const frame);

    // Make the velocity on the grid to be divergence free
    void pressureProjection(bool print);
    void pressureProjection4(bool print);

    void updateEmitter();

    void createDirichletVelocity();
    void createDirichletVelocity2();
    void createDirichletVelocity4();
    void createVCurr4(bool print);
    void createDensityCurr4();
    void createFlags4(bool print);
    void createInteriorPressure4();

    void applyDirichletVelocity(Vec3SGrid& vecGrid, int frame);
    // going over mDirichletVelocity
    void applyDirichletVelocity4(Vec3SGrid& vecGrid, int frame);

    void swap();

    void addGravity(float const dt);

    void advectDensity(float const dt);

    void advectVelocity(float const dt, int const frame);

    void writeVDBs(int const frame);

    void writeVDBsDebug(int const frame);

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
                else if (neighbor.x() - 1 == ijk.x() /* right x-face */) {
                    diagonal -= 1.0;
                    source -= dirichletBC;
                }
                else if (neighbor.y() + 1 == ijk.y() /* bottom y-face */) {
                    diagonal -= 1.0;
                    source -= dirichletBC;
                }
                else if (neighbor.y() - 1 == ijk.y() /* top y-face */) {
                    diagonal -= 1.0;
                    source -= dirichletBC;
                }
                else if (neighbor.z() + 1 == ijk.z() /* back z-face */) {
                    diagonal -= 1.0;
                    source -= dirichletBC;
                }
                else if (neighbor.z() - 1 == ijk.z() /* front z-face */) {
                    diagonal -= 1.0;
                    source -= dirichletBC;
                }
            }
        }

        Vec3SGrid::ConstPtr dirichletVelocity;
        FloatGrid::ConstPtr dirichletPressure;
        float voxelSize;
    };

    struct BoundaryOp4 {
        BoundaryOp4(Int32Grid::ConstPtr flags,
                    Vec3SGrid::ConstPtr dirichletVelocity,
                    FloatGrid::ConstPtr dirichletPressure,
                    float const voxelSize) :
                    flags(flags),
                    dirichletVelocity(dirichletVelocity),
                    dirichletPressure(dirichletPressure),
                    voxelSize(voxelSize) {}

        void operator()(const openvdb::Coord& ijk,
                        const openvdb::Coord& neighbor,
                        double& source,
                        double& diagonal) const
        {
            float const dirichletBC = 0.f;
            int flag = flags->tree().getValue(neighbor);
            bool isNeumannPressure = (flag == 0);
            bool isDirichletPressure = (flag == 4);
            auto vNgbr = Vec3s::zero(); //dirichletVelocity->tree().getValue(neighbor);

            // TODO: Double check this:
            if (isNeumannPressure) {
                double delta = 0.0;
                // Neumann pressure from bbox
                if (neighbor.x() + 1 == ijk.x() /* left x-face */) {
                    delta += vNgbr[0];

                    // std::cout << "ijk = " << ijk << "\tneighbor = " << neighbor
                    //           << "\tneumann = " << isNeumannPressure
                    //           << "\tvngbr = " << vNgbr
                    //           << "\tleft x face"
                    //           << "\tdelta = " << delta
                    //           << std::endl;
                }
                if (neighbor.x() - 1 == ijk.x() /* right x-face */) {
                    delta -= vNgbr[0];
                    // std::cout << "ijk = " << ijk << "\tneighbor = " << neighbor
                    //           << "\tneumann = " << isNeumannPressure
                    //           << "\tvngbr = " << vNgbr
                    //           << "\tright x face"
                    //           << "\tdelta = " << delta
                    //           << std::endl;
                }
                if (neighbor.y() + 1 == ijk.y() /* bottom y-face */) {
                    delta += vNgbr[1];
                    // std::cout << "ijk = " << ijk << "\tneighbor = " << neighbor
                    //           << "\tneumann = " << isNeumannPressure
                    //           << "\tvngbr = " << vNgbr
                    //           << "\tbottom y face"
                    //           << "\tdelta = " << delta
                    //           << std::endl;
                }
                if (neighbor.y() - 1 == ijk.y() /* top y-face */) {
                    delta -= vNgbr[1];
                    // std::cout << "ijk = " << ijk << "\tneighbor = " << neighbor
                    //           << "\tneumann = " << isNeumannPressure
                    //           << "\tvngbr = " << vNgbr
                    //           << "\ttop y face"
                    //           << "\tdelta = " << delta
                    //           << std::endl;
                }
                if (neighbor.z() + 1 == ijk.z() /* back z-face */) {
                    delta +=  vNgbr[2];
                    // std::cout << "ijk = " << ijk << "\tneighbor = " << neighbor
                    //           << "\tneumann = " << isNeumannPressure
                    //           << "\tvngbr = " << vNgbr
                    //           << "\tback z face"
                    //           << "\tdelta = " << delta
                    //           << std::endl;
                }
                if (neighbor.z() - 1 == ijk.z() /* front z-face */) {
                    delta -=  vNgbr[2];
                    // std::cout << "ijk = " << ijk << "\tneighbor = " << neighbor
                    //           << "\tneumann = " << isNeumannPressure
                    //           << "\tvngbr = " << vNgbr
                    //           << "\tfront z face"
                    //           << "\tdelta = " << delta
                    //           << std::endl;
                }
                // Note: in the SOP_OpenVDB_Remove_Divergence, we need to multiply
                // this by 0.5, because the gradient that's used is using
                // central-differences in a collocated grid, instead of the staggered one.
                source += delta / voxelSize;
            } else if (isDirichletPressure) {
                std::cout << "dirichlet" << std::endl;
                // Dirichlet pressure
                if (neighbor.x() + 1 == ijk.x() /* left x-face */) {
                    diagonal -= 1.0;
                    source -= dirichletBC;
                }
                else if (neighbor.x() - 1 == ijk.x() /* right x-face */) {
                    diagonal -= 1.0;
                    source -= dirichletBC;
                }
                else if (neighbor.y() + 1 == ijk.y() /* bottom y-face */) {
                    diagonal -= 1.0;
                    source -= dirichletBC;
                }
                else if (neighbor.y() - 1 == ijk.y() /* top y-face */) {
                    diagonal -= 1.0;
                    source -= dirichletBC;
                }
                else if (neighbor.z() + 1 == ijk.z() /* back z-face */) {
                    diagonal -= 1.0;
                    source -= dirichletBC;
                }
                else if (neighbor.z() - 1 == ijk.z() /* front z-face */) {
                    diagonal -= 1.0;
                    source -= dirichletBC;
                }
            }
        }

        Int32Grid::ConstPtr flags;
        Vec3SGrid::ConstPtr dirichletVelocity;
        FloatGrid::ConstPtr dirichletPressure;
        float voxelSize;
    };



    // Apply Gravity Functor. Meant to be used with
    // foreach in LeafManager
    struct ApplyGravityOp
    {
        ApplyGravityOp(float const dt, Vec3s const gravity, Int32Grid::Ptr flags) :
            dt(dt), gravity(gravity), flags(flags) {}

        template <typename T>
        void operator()(T &leaf, size_t) const
        {
            auto flagsAcc = flags->getConstAccessor(); 
            for (typename T::ValueOnIter iter = leaf.beginValueOn(); iter; ++iter) {
                auto ijk = iter.getCoord();
                Vec3s newVal = *iter  + dt * gravity;
                //if (flagsAcc.getValue(ijk) == 1) {
                    iter.setValue(newVal);
                //}
            }
        }

        Vec3s const gravity;
        float const dt;
        Int32Grid::Ptr flags;
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
    Vec3s mGravity = Vec3s(0.f, -48.f, 0.f);
    Vec3s mPushVelocity = Vec3s(1.f, 0.f, 0.f);
    math::Transform::Ptr mXform;
    float mPadding;
    Coord mMin;
    Coord mMax;
    Coord mMaxStaggered;

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
    BoolGrid::Ptr mCollocatedMaskGrid;
    Int32Grid::Ptr mFlags;
    FloatGrid::Ptr mSphere;
    BoolGrid::Ptr mInteriorPressure;
};


SmokeSolver::SmokeSolver(float const voxelSize) : mVoxelSize(voxelSize)
{
    initialize4();
}

 void
 SmokeSolver::createFlags4(bool print)
 {
    mFlags = Int32Grid::create(/* bg = */ 0); // Neumann pressure
    mFlags->denseFill(CoordBBox(mMin, mMax), /* value = */ 1, /* active = */ true);
    mFlags->setTransform(mXform);
    mFlags->setName("flags");
    auto flagsAcc = mFlags->getAccessor();
    auto sphereAcc = mSphere->getAccessor();
    for (auto iter = mFlags->beginValueOn(); iter; ++iter) {
        math::Coord ijk = iter.getCoord();

        if (/* ijk[0] == mMax[0] */ ijk[1] == mMax[1]) {
            flagsAcc.setValue(ijk, 4); // Dirichlet
        }
        if (ijk[0] == mMin[0] /* left face */ ||
            ijk[1] == mMin[1] /* bottom face */ ||
            /*ijk[1] == mMax[1]*/ /* top face */ ijk[0] == mMax[0]||
            ijk[2] == mMin[2] /* back face */ ||
            ijk[2] == mMax[2] /* front face */ ||
            sphereAcc.getValue(ijk) < 0) {
            flagsAcc.setValue(ijk, 0); // Neumann
        }
        if (print) {
            std::cout << "flags" << ijk << " = " << flagsAcc.getValue(ijk) << std::endl;
        }
    }
    std::ostringstream ostr;
    ostr << "flags.vdb";
    std::cerr << "\tWriting " << ostr.str() << std::endl;
    openvdb::io::File file(ostr.str());
    openvdb::GridPtrVec grids;
    grids.push_back(mFlags);
    file.write(grids);
    file.close();
 }


void
SmokeSolver::createInteriorPressure4()
{
    mInteriorPressure = BoolGrid::create(false);
    mInteriorPressure->denseFill(CoordBBox(mMin, mMax), /* value = */ true, /* active = */ true);
    mInteriorPressure->setTransform(mXform);
    mInteriorPressure->setName("interior_pressure");

    auto flagsAcc = mFlags->getConstAccessor();
    for (auto iter = mInteriorPressure->beginValueOn(); iter; ++iter) {
        math::Coord ijk = iter.getCoord();
        if (flagsAcc.getValue(ijk) != 1) {
            iter.setValueOff();
        }
    }

    // std::cout << "\ninterior pressure" << std::endl;
    // auto ipAcc = mInteriorPressure->getAccessor();
    // for (auto iter = mInteriorPressure->beginValueOn(); iter; ++iter) {
    //     std::cout << "int pres ijk = " << iter.getCoord() << std::endl;
    // }
}

 void
 SmokeSolver::createDensityCurr4()
 {
    mDensityCurr = FloatGrid::create(/*bg = */0.f);
    mDensityCurr->setTransform(mXform);
    mDensityCurr->setName("density_curr");
    mDensityCurr->denseFill(CoordBBox(mMin, mMax), /* value = */ 0.f, /* active = */ true);
    mDensityCurr->topologyUnion(*mFlags);
    mDensityCurr->topologyIntersection(*mFlags);
 }


 void
 SmokeSolver::createVCurr4(bool print)
 {
    mVCurr = Vec3SGrid::create(/* bg = */ Vec3s::zero()); // Neumann pressure
    mVCurr->setGridClass(GRID_STAGGERED);
    // mVCurr->denseFill(CoordBBox(mMin, mMaxStaggered), /* value = */ Vec3s(1.0, 0.f, 0.f), /* active = */ true);
    // 5.15 am
    mVCurr->denseFill(CoordBBox(mMin, mMax), /* value = */ mPushVelocity, /* active = */ true);
    mVCurr->setTransform(mXform);
    mVCurr->setName("vel_curr");

    if (print) {
        std::cout << "\ncreate vcurr4 velocity" << std::endl;
        auto velAcc = mVCurr->getAccessor();
        for (auto iter = mVCurr->beginValueOn(); iter; ++iter) {
            auto ijk = iter.getCoord();
            std::cout << "vel" << ijk  << " = " << velAcc.getValue(ijk) << std::endl;
        }
    }

 }


 void
 SmokeSolver::initialize4()
 {
    std::cout << "initialize 4" << std::endl;
    using BBox = math::BBox<Vec3s>;
    mPushVelocity = Vec3s(1.f, 0.f, 0.f);
    mXform = math::Transform::createLinearTransform(mVoxelSize);
    float const padding = 1.f * mVoxelSize;
    mPadding = padding;
    float const centerY = 3.f;
    auto minBBox = Vec3s(0.f, 0.f, 0.f);
    // works:
    auto maxBBox = Vec3s(3 * mVoxelSize, 3 * mVoxelSize, 3 * mVoxelSize);
    //auto maxBBox = Vec3s(1.f + mVoxelSize, 0.4f + mVoxelSize, 0.4f + mVoxelSize);
    // auto maxBBox = Vec3s(7.f + mVoxelSize, 3.f + mVoxelSize, 3.f + mVoxelSize);
    Coord minBBoxIntrCoord = mXform->worldToIndexNodeCentered(minBBox);
    Coord maxBBoxIntrCoord = mXform->worldToIndexNodeCentered(maxBBox);
    mMin = minBBoxIntrCoord;
    mMax = maxBBoxIntrCoord;
    mMaxStaggered = mMax + Coord(1);

    const float radius = 0.3 * maxBBox[1];
    const openvdb::Vec3f center(100 + 2.5f / 7.f * maxBBox[0], 0.5f * maxBBox[1], 0.5f * maxBBox[2]);
    mSphere = tools::createLevelSetSphere<openvdb::FloatGrid>(radius, center, mVoxelSize, 2 /* width */);

    std::cout << "mMax = " << mMax << "\tmMaxStaggered = " << mMaxStaggered << std::endl;
    mCollocatedMaskGrid = BoolGrid::create(false /* background */);
    mCollocatedMaskGrid->denseFill(CoordBBox(minBBoxIntrCoord, maxBBoxIntrCoord), /* value = */ true, /* active = */ true);
    mCollocatedMaskGrid->setTransform(mXform);
    mCollocatedMaskGrid->setName("domain_mask");

    // Create an emitter and an emitter velocity
    auto minEmtW = Vec3s(0.f, 2.5f, 2.5f);
    auto maxEmtW = Vec3s(2 * mVoxelSize, 3.5f, 3.5f);
    Coord minEmtCoord = mXform->worldToIndexNodeCentered(minEmtW);
    Coord maxEmtCoord = mXform->worldToIndexNodeCentered(maxEmtW);
    mEmitter = FloatGrid::create(/*bg = */0.f);
    mEmitter->denseFill(CoordBBox(minEmtCoord, maxEmtCoord), /* value = */ 2.0, /*active = */ true);
    mEmitter->setTransform(mXform);
    mEmitter->setName("emitter");

    createFlags4(false);
    createInteriorPressure4();
    createVCurr4(false);
    // createDensityCurr4();
    createDirichletVelocity4();

    // updateEmitter();
    applyDirichletVelocity4(*mVCurr, -1);

    addGravity(1.f/24.f);
    applyDirichletVelocity4(*mVCurr, -1);
    pressureProjection4(false);

    writeVDBsDebug(-1);
    exit(0);
 }

 void
 SmokeSolver::pressureProjection4(bool print)
 {
    std::cout << "pressure projection 4" << std::endl;
    using TreeType = FloatTree;
    using ValueType = TreeType::ValueType;
    using MaskGridType = BoolGrid;
    using PCT = openvdb::math::pcg::JacobiPreconditioner<openvdb::tools::poisson::LaplacianMatrix>;


    ValueType const zero = zeroVal<ValueType>();
    double const epsilon = math::Delta<ValueType>::value();

    mDivBefore = tools::divergence(*mVCurr);
    mDivBefore->topologyIntersection(*mInteriorPressure);
    mDivBefore->setName("div_before");

    float divBefore = 0.f;
    auto divBeforeAcc = mDivBefore->getAccessor();
    auto flagAcc = mFlags->getAccessor();
    auto vCurrAcc = mVCurr->getAccessor();
    //for (auto iter = mDivBefore->beginValueOn(); iter; ++iter) {
    for (int kk = 1; kk <= 2; ++kk)
    for (int jj = 1; jj <= 2; ++jj)
    for (int ii = 1; ii <= 2; ++ii) {
        math::Coord ijk(ii, jj, kk); //= iter.getCoord();
        auto ip1jk = ijk.offsetBy(1, 0, 0);
        auto ijp1k = ijk.offsetBy(0, 1, 0);
        auto ijkp1 = ijk.offsetBy(0, 0, 1);
        auto val = divBeforeAcc.getValue(ijk);

        Vec3s vdown(vCurrAcc.getValue(ijk));
        Vec3s vup(vCurrAcc.getValue(ip1jk)[0],
                  vCurrAcc.getValue(ijp1k)[1],
                  vCurrAcc.getValue(ijkp1)[2]);
                  if (print) {
                    std::cout << "div before " << ijk << " = " << val
                              << " vdown = " << vdown
                              << " vup = " << vup
                              << "flags[ijk]" << flagAcc.getValue(ijk)
                              << std::endl;
                  }

        if (std::abs(val) > std::abs(divBefore)) {
            divBefore = val;
        }
    }

    std::cout << "\t== divergence before pp = " << divBefore << std::endl;

    math::pcg::State state = math::pcg::terminationDefaults<ValueType>();
    state.iterations = 100000;
    state.relativeError = state.absoluteError = epsilon;
    SmokeSolver::BoundaryOp4 bop(mFlags, mVCurr, mDirichletPressure, mVoxelSize);
    util::NullInterrupter interrupter;
    FloatTree::Ptr fluidPressure = tools::poisson::solveWithBoundaryConditionsAndPreconditioner<PCT>(
        mDivBefore->tree(), mInteriorPressure->tree(), bop, state, interrupter, /*staggered=*/true);


    std::cout << "Projection Success: " << state.success << "\n";
    std::cout << "Iterations: " << state.iterations << "\n";
    std::cout << "Relative error: " << state.relativeError << "\n";
    std::cout << "Absolute error: " << state.absoluteError << "\n";

    // Note: need to dilate in order to do one-sided difference
    // because we use a staggered grid velocity field.
    FloatGrid::Ptr fluidPressureGrid = FloatGrid::create(fluidPressure);
    // tools::dilateActiveValues(*fluidPressure, /*iterations=*/1, tools::NN_FACE, tools::IGNORE_TILES);

    fluidPressureGrid->setTransform(mXform);
    mPressure = fluidPressureGrid->copy();
    mPressure->setName("pressure");

    auto pressureAcc = fluidPressureGrid->getConstAccessor();
    auto flagsAcc = mFlags->getConstAccessor();


    for (auto iter = mPressure->beginValueOn(); iter; ++iter) {
        auto ijk = iter.getCoord();
        std::cout << "p" << ijk << " = " << pressureAcc.getValue(ijk) << std::endl;
    }



    // Note: I'm modifying vCurr
    for (auto iter = mVCurr->beginValueOn(); iter; ++iter) {
        auto ijk = iter.getCoord();
        auto im1jk = ijk.offsetBy(-1, 0, 0);
        auto ijm1k = ijk.offsetBy(0, -1, 0);
        auto ijkm1 = ijk.offsetBy(0, 0, -1);

        // Only updates velocity if it is a face of fluid cell

        if (flagsAcc.getValue(ijk) == 1 ||
            flagsAcc.getValue(im1jk) == 1 || 
            flagsAcc.getValue(ijm1k) == 1 || 
            flagsAcc.getValue(ijkm1) == 1) {
            Vec3s gradijk;
            gradijk[0] = pressureAcc.getValue(ijk) - pressureAcc.getValue(ijk.offsetBy(-1, 0, 0));
            gradijk[1] = pressureAcc.getValue(ijk) - pressureAcc.getValue(ijk.offsetBy(0, -1, 0));
            gradijk[2] = pressureAcc.getValue(ijk) - pressureAcc.getValue(ijk.offsetBy(0, 0, -1));
            auto val = vCurrAcc.getValue(ijk) - gradijk * mVoxelSize;
            vCurrAcc.setValue(ijk, val);
        }
    }

    applyDirichletVelocity4(*mVCurr, -2);
    mDivAfter = tools::divergence(*mVCurr);
    mDivAfter->setName("div_after");
    (mDivAfter->tree()).topologyIntersection(mInteriorPressure->tree());
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

    // std::ostringstream ostr;
    // ostr << "debug_divergence.vdb";
    // std::cerr << "\tWriting " << ostr.str() << std::endl;
    // openvdb::io::File file(ostr.str());
    // openvdb::GridPtrVec grids;
    // grids.push_back(mDivBefore);
    // grids.push_back(mDivAfter);
    // grids.push_back(mPressure);
    // grids.push_back(mVCurr);
    // grids.push_back(mDensityCurr);
    // file.write(grids);
    // file.close();
 }

 void
 SmokeSolver::createDirichletVelocity4() {
    std::cout << "create dirichlet velocity 4" << std::endl;

    mDirichletVelocity = Vec3SGrid::create(/* bg = */ Vec3s(0.f, 0.f, 0.f));
    mDirichletVelocity->setName("dirichlet_velocity");
    mDirichletVelocity->setTransform(mXform);
    mDirichletVelocity->topologyUnion(*mVCurr);
    mDirichletVelocity->setGridClass(GRID_STAGGERED);

    auto flagsAcc = mFlags->getConstAccessor();
    auto drcAcc = mDirichletVelocity->getAccessor();

    for (auto iter = mDirichletVelocity->beginValueOn(); iter; ++iter) {
        auto ijk = iter.getCoord();
        auto im1jk = ijk.offsetBy(-1, 0, 0);
        auto ijm1k = ijk.offsetBy(0, -1, 0);
        auto ijkm1 = ijk.offsetBy(0, 0, -1);

        Vec3s val;

        if (flagsAcc.getValue(ijk) == 0) {
            // is a full Neumann
            val = mPushVelocity;
        } else {
            if(flagsAcc.getValue(im1jk) == 0) {
                // neighboring a Neumann pressure in the x face
                val[0] = mPushVelocity[0];
            }
            if(flagsAcc.getValue(ijm1k) == 0) {
                // neighboring a Neumann pressure in the y face
                val[1] = mPushVelocity[1];
            }
            if(flagsAcc.getValue(ijkm1) == 0) {
                // neighboring a Neumann pressure in the z face
                val[2] = mPushVelocity[2];
            }
        }
        drcAcc.setValue(ijk, val);
    }

    // mDirichletVelocity = Vec3SGrid::create(/* bg = */ Vec3s(0.f, 0.f, 0.f));
    // mDirichletVelocity->setName("dirichlet_velocity");
    // mDirichletVelocity->setTransform(mXform);
    // mDirichletVelocity->topologyUnion(*mFlags);

    // Vec3s const pushVelocity(0.f, 0.f, 0.f);

    // auto flagsAcc = mFlags->getConstAccessor();
    // auto drcAcc = mDirichletVelocity->getAccessor();

    // for (auto iter = mDirichletVelocity->beginValueOn(); iter; ++iter) {
    //     auto ijk = iter.getCoord();
    //     // Neumann voxel
    //     if (flagsAcc.getValue(ijk) == 0) {
    //         if (ijk[0] == 0) {
    //             drcAcc.setValue(ijk, pushVelocity);
    //         } else {
    //             drcAcc.setValue(ijk, Vec3s::zero());
    //         }
    //     } else {
    //         iter.setValueOff();
    //     }
    // }
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
SmokeSolver::applyDirichletVelocity4(Vec3SGrid& vecGrid, int frame)
{
    std::cout << "apply dirichlet velocity begins" << std::endl;
    auto vecAcc = vecGrid.getAccessor();
    auto drcAcc = mDirichletVelocity->getConstAccessor();
    auto flagsAcc = mFlags->getConstAccessor();
    for (auto iter = vecGrid.beginValueOn(); iter; ++iter) {
        auto ijk = iter.getCoord();
        auto im1jk = ijk.offsetBy(-1, 0, 0);
        auto ijm1k = ijk.offsetBy(0, -1, 0);
        auto ijkm1 = ijk.offsetBy(0, 0, -1);

        Vec3s val = *iter;

        if (flagsAcc.getValue(ijk) == 0) {
            // is a full Neumann
            val = drcAcc.getValue(ijk);
        } else {
            if(flagsAcc.getValue(im1jk) == 0) {
                // neighboring a Neumann pressure in the x face
                // not looking at getValue(im1jk),
                // because that's how we set dirichlet velocity being staggered
                val[0] = drcAcc.getValue(ijk)[0];
            }
            if(flagsAcc.getValue(ijm1k) == 0) {
                // neighboring a Neumann pressure in the y face
                val[1] = drcAcc.getValue(ijk)[1];
            }
            if(flagsAcc.getValue(ijkm1) == 0) {
                // neighboring a Neumann pressure in the z face
                val[2] = drcAcc.getValue(ijk)[2];
            }
        }
        iter.setValue(val);
    }
    bool print = false;
    if (print) {
    std::cout << "\nafter boundary correction" << std::endl;
    auto velAcc = mVCurr->getAccessor();
    for (auto iter = mVCurr->beginValueOn(); iter; ++iter) {
        auto ijk = iter.getCoord();
        std::cout << "bdry corrected vel" << ijk  << " = " << velAcc.getValue(ijk) << std::endl;
    }
    }
    

    // // First attempt
    // auto vecAcc = vecGrid.getAccessor();
    // auto drcAcc = mDirichletVelocity->getConstAccessor();
    // bool print = true;
    // for (auto iter = mDirichletVelocity->beginValueOn(); iter; ++iter) {
    //     auto ijk = iter.getCoord();
    //     auto ip1jk = ijk.offsetBy(1, 0, 0);
    //     auto ijp1k = ijk.offsetBy(0, 1, 0);
    //     auto ijkp1 = ijk.offsetBy(0, 0, 1);

    //     Vec3s val = *iter;

    //     vecAcc.setValue(ijk, *iter);
    //     if (vecAcc.isValueOn(ip1jk)) {
    //         auto oldVel = vecAcc.getValue(ip1jk);
    //         auto newVel = Vec3s(val[0], oldVel[1], oldVel[2]);
    //         vecAcc.setValue(ip1jk, newVel);
    //     }
    //     if (vecAcc.isValueOn(ijp1k)) {
    //         auto oldVel = vecAcc.getValue(ijp1k);
    //         auto newVel = Vec3s(oldVel[0], val[1], oldVel[2]);
    //         vecAcc.setValue(ijp1k, newVel);
    //     }
    //     if (vecAcc.isValueOn(ijkp1)) {
    //         auto oldVel = vecAcc.getValue(ijkp1);
    //         auto newVel = Vec3s(oldVel[0], oldVel[1], val[2]);
    //         vecAcc.setValue(ijkp1, newVel);
    //     }
    // }


    // std::cout << "apply dirichlet velocity begins" << std::endl;
    // auto vecAcc = vecGrid.getAccessor();
    // auto drcAcc = mDirichletVelocity->getConstAccessor();
    // auto flagsAcc = mFlags->getConstAccessor();
    // bool print = true;
    // for (auto iter = vecGrid.beginValueOn(); iter; ++iter) {
    //     auto ijk = iter.getCoord();
    //     auto im1jk = ijk.offsetBy(-1, 0, 0);
    //     auto ijm1k = ijk.offsetBy(0, -1, 0);
    //     auto ijkm1 = ijk.offsetBy(0, 0, -1);
    //     auto drcVel = drcAcc.getValue(ijk);

    //     if (flagsAcc.getValue(ijk) == 0) {
    //         iter.setValue(drcVel);
    //     }

    //     if (flagsAcc.getValue(im1jk) == 0) {
    //         auto val = drcAcc.getValue(im1jk);
    //         auto oldVel = vecAcc.getValue(im1jk);
    //         auto newVel = Vec3s(val[0], oldVel[1], oldVel[2]);
    //         vecAcc.setValue(im1jk, newVel);
    //     }

    //     if (flagsAcc.getValue(ijm1k) == 0) {
    //         auto val = drcAcc.getValue(ijm1k);
    //         auto oldVel = vecAcc.getValue(ijm1k);
    //         auto newVel = Vec3s(oldVel[0], val[1], oldVel[2]);
    //         vecAcc.setValue(ijm1k, newVel);
    //     }


    //     if (flagsAcc.getValue(ijkm1) == 0) {
    //         auto val = drcAcc.getValue(ijkm1);
    //         auto oldVel = vecAcc.getValue(ijkm1);
    //         auto newVel = Vec3s(oldVel[0], oldVel[1], val[2]);
    //         vecAcc.setValue(ijkm1, newVel);

    //     }
    // }
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
SmokeSolver::createDirichletVelocity()
{
    auto minDVLft = Vec3s(-mPadding, -mPadding, -mPadding);
    auto maxDVLft = Vec3s(2 * mVoxelSize, mMax[1], mMax[2]);
    Coord minDVLftCoord = mXform->worldToIndexNodeCentered(minDVLft);
    Coord maxDVLftCoord = mXform->worldToIndexNodeCentered(maxDVLft);
    mDirichletVelocity = Vec3SGrid::create(/* bg = */ Vec3s(0.f, 0.f, 0.f));
    mDirichletVelocity->denseFill(CoordBBox(minDVLftCoord, maxDVLftCoord), /* value = */ Vec3s(1.f, 0.f, 0.f), /*active = */ true);
    mDirichletVelocity->setGridClass(GRID_STAGGERED);
    mDirichletVelocity->setTransform(mXform);
    mDirichletVelocity->setName("dirichlet_velocity");

    auto minDVBck = Vec3s(2 * mVoxelSize, 0.f, -mPadding);
    auto maxDVBck = Vec3s(14.f + mPadding, 6.f, 0.f);
    Coord minDVBckCoord = mXform->worldToIndexNodeCentered(minDVBck);
    Coord maxDVBckCoord = mXform->worldToIndexNodeCentered(maxDVBck);
    Vec3SGrid::Ptr bck = Vec3SGrid::create(/* bg = */ Vec3s(0.f, 0.f, 0.f));
    bck->setGridClass(GRID_STAGGERED);
    bck->denseFill(CoordBBox(minDVBckCoord, maxDVBckCoord), /* value = */ Vec3s(0.f, 0.f, 0.f), /*active = */ true);

    auto minDVFrt = Vec3s(2 * mVoxelSize, 0.f, 6.f);
    auto maxDVFrt = Vec3s(14.f + mPadding, 6.f, 6.f + mPadding);
    Coord minDVFrtCoord = mXform->worldToIndexNodeCentered(minDVFrt);
    Coord maxDVFrtCoord = mXform->worldToIndexNodeCentered(maxDVFrt);
    Vec3SGrid::Ptr frt = Vec3SGrid::create(/* bg = */ Vec3s(0.f, 0.f, 0.f));
    frt->setGridClass(GRID_STAGGERED);
    frt->denseFill(CoordBBox(minDVFrtCoord, maxDVFrtCoord), /* value = */ Vec3s(0.f, 0.f, 0.f), /*active = */ true);

    auto minDVTop = Vec3s(2 * mVoxelSize, 6.f, 0.f);
    auto maxDVTop = Vec3s(14.f + mPadding, 6.f + mPadding, 6.f);
    Coord minDVTopCoord = mXform->worldToIndexNodeCentered(minDVTop);
    Coord maxDVTopCoord = mXform->worldToIndexNodeCentered(maxDVTop);
    Vec3SGrid::Ptr top = Vec3SGrid::create(/* bg = */ Vec3s(0.f, 0.f, 0.f));
    top->setGridClass(GRID_STAGGERED);
    top->denseFill(CoordBBox(minDVTopCoord, maxDVTopCoord), /* value = */ Vec3s(0.f, 0.f, 0.f), /*active = */ true);

    auto minDVBtm = Vec3s(2 * mVoxelSize, -mPadding, 0.f);
    auto maxDVBtm = Vec3s(14.f + mPadding, 0.f, 6.f);
    Coord minDVBtmCoord = mXform->worldToIndexNodeCentered(minDVBtm);
    Coord maxDVBtmCoord = mXform->worldToIndexNodeCentered(maxDVBtm);
    Vec3SGrid::Ptr btm = Vec3SGrid::create(/* bg = */ Vec3s(0.f, 0.f, 0.f));
    btm->setGridClass(GRID_STAGGERED);
    btm->denseFill(CoordBBox(minDVBtmCoord, maxDVBtmCoord), /* value = */ Vec3s(0.f, 0.f, 0.f), /*active = */ true);

    mDirichletVelocity->topologyUnion(*bck);
    mDirichletVelocity->topologyUnion(*frt);
    mDirichletVelocity->topologyUnion(*top);
    mDirichletVelocity->topologyUnion(*btm);
}



void
SmokeSolver::createDirichletVelocity2()
{
    BoolGrid::Ptr interiorMaskGrid = BoolGrid::create(false);
    interiorMaskGrid->topologyUnion(*mCollocatedMaskGrid);
    tools::erodeActiveValues(interiorMaskGrid->tree(), /*iterations=*/1, tools::NN_FACE, tools::IGNORE_TILES);
    mCollocatedMaskGrid->setName("interior_mask_grid");

    mDirichletVelocity = Vec3SGrid::create(/* bg = */ Vec3s(0.f, 0.f, 0.f));
    mDirichletVelocity->setName("dirichlet_velocity");
    mDirichletVelocity->setTransform(mXform);
    mDirichletVelocity->setGridClass(GRID_STAGGERED);
    mDirichletVelocity->denseFill(CoordBBox(mMin, Coord(1, mMax[1], mMax[2])), /* value = */ Vec3s(1.f, 0.f, 0.f), /*active = */ true);
    mDirichletVelocity->topologyUnion(*mCollocatedMaskGrid);
    tools::dilateActiveValues(mDirichletVelocity->tree(), /*iterations=*/1, tools::NN_FACE, tools::IGNORE_TILES);
    mDirichletVelocity->topologyDifference(*interiorMaskGrid);

    auto dirichletAcc = mDirichletVelocity->getAccessor();
    for (auto iter = mDirichletVelocity->beginValueOn(); iter; ++iter) {
        math::Coord ijk = iter.getCoord();
        if (ijk[0] >= mMax[0])
        dirichletAcc.setValueOff(ijk);
    }
    std::ostringstream ostr;
    ostr << "debug_dirichlet_velocity" << ".vdb";
    std::cerr << "Writing " << ostr.str() << std::endl;
    openvdb::io::File file(ostr.str());
    openvdb::GridPtrVec grids;
    grids.push_back(mDirichletVelocity);
    grids.push_back(mCollocatedMaskGrid);
    file.write(grids);
    file.close();

    //mDirichletVelocity->denseFill(CoordBBox(minDVLftCoord, maxDVLftCoord), /* value = */ Vec3s(1.f, 0.f, 0.f), /*active = */ true);
    //mDirichletVelocity->setGridClass(GRID_STAGGERED);
    //mDirichletVelocity->setTransform(mXform);
    //mDirichletVelocity->setName("dirichlet_velocity");

    //auto minDVBck = Vec3s(2 * mVoxelSize, 0.f, -mPadding);
    //auto maxDVBck = Vec3s(14.f + mPadding, 6.f, 0.f);
    //Coord minDVBckCoord = mXform->worldToIndexNodeCentered(minDVBck);
    //Coord maxDVBckCoord = mXform->worldToIndexNodeCentered(maxDVBck);
    //Vec3SGrid::Ptr bck = Vec3SGrid::create(/* bg = */ Vec3s(0.f, 0.f, 0.f));
    //bck->setGridClass(GRID_STAGGERED);
    //bck->denseFill(CoordBBox(minDVBckCoord, maxDVBckCoord), /* value = */ Vec3s(0.f, 0.f, 0.f), /*active = */ true);

    //auto minDVFrt = Vec3s(2 * mVoxelSize, 0.f, 6.f);
    //auto maxDVFrt = Vec3s(14.f + mPadding, 6.f, 6.f + mPadding);
    //Coord minDVFrtCoord = mXform->worldToIndexNodeCentered(minDVFrt);
    //Coord maxDVFrtCoord = mXform->worldToIndexNodeCentered(maxDVFrt);
    //Vec3SGrid::Ptr frt = Vec3SGrid::create(/* bg = */ Vec3s(0.f, 0.f, 0.f));
    //frt->setGridClass(GRID_STAGGERED);
    //frt->denseFill(CoordBBox(minDVFrtCoord, maxDVFrtCoord), /* value = */ Vec3s(0.f, 0.f, 0.f), /*active = */ true);

    //auto minDVTop = Vec3s(2 * mVoxelSize, 6.f, 0.f);
    //auto maxDVTop = Vec3s(14.f + mPadding, 6.f + mPadding, 6.f);
    //Coord minDVTopCoord = mXform->worldToIndexNodeCentered(minDVTop);
    //Coord maxDVTopCoord = mXform->worldToIndexNodeCentered(maxDVTop);
    //Vec3SGrid::Ptr top = Vec3SGrid::create(/* bg = */ Vec3s(0.f, 0.f, 0.f));
    //top->setGridClass(GRID_STAGGERED);
    //top->denseFill(CoordBBox(minDVTopCoord, maxDVTopCoord), /* value = */ Vec3s(0.f, 0.f, 0.f), /*active = */ true);

    //auto minDVBtm = Vec3s(2 * mVoxelSize, -mPadding, 0.f);
    //auto maxDVBtm = Vec3s(14.f + mPadding, 0.f, 6.f);
    //Coord minDVBtmCoord = mXform->worldToIndexNodeCentered(minDVBtm);
    //Coord maxDVBtmCoord = mXform->worldToIndexNodeCentered(maxDVBtm);
    //Vec3SGrid::Ptr btm = Vec3SGrid::create(/* bg = */ Vec3s(0.f, 0.f, 0.f));
    //btm->setGridClass(GRID_STAGGERED);
    //btm->denseFill(CoordBBox(minDVBtmCoord, maxDVBtmCoord), /* value = */ Vec3s(0.f, 0.f, 0.f), /*active = */ true);

    //mDirichletVelocity->topologyUnion(*bck);
    //mDirichletVelocity->topologyUnion(*frt);
    //mDirichletVelocity->topologyUnion(*top);
    //mDirichletVelocity->topologyUnion(*btm);
}


void
SmokeSolver::initialize() {
    std::cout << "initialize " << std::endl;
    using BBox = math::BBox<Vec3s>;

    mXform = math::Transform::createLinearTransform(mVoxelSize);
    float const padding = 4.f * mVoxelSize;
    mPadding = padding;
    float const centerY = 3.f;
    auto minBBox = Vec3s(0.f, 0.f, 0.f);
    auto maxBBox = Vec3s(14.f, 6.f, 6.f);
    Coord minBBoxIntrCoord = mXform->worldToIndexNodeCentered(minBBox);
    Coord maxBBoxIntrCoord = mXform->worldToIndexNodeCentered(maxBBox);
    mMin = minBBoxIntrCoord;
    mMax = maxBBoxIntrCoord;
    mMaxStaggered = mMax + Coord(1);
    mCollocatedMaskGrid = BoolGrid::create(false /* background */);
    mCollocatedMaskGrid->denseFill(CoordBBox(minBBoxIntrCoord, maxBBoxIntrCoord), /* value = */ true, /* active = */ true);
    mCollocatedMaskGrid->setTransform(mXform);
    mCollocatedMaskGrid->setName("domain_mask");

    // Create an emitter and an emitter velocity
    auto minEmtW = Vec3s(0.f, 2.5f, 2.5f);
    auto maxEmtW = Vec3s(2 * mVoxelSize, 3.5f, 3.5f);
    Coord minEmtCoord = mXform->worldToIndexNodeCentered(minEmtW);
    Coord maxEmtCoord = mXform->worldToIndexNodeCentered(maxEmtW);
    mEmitter = FloatGrid::create(/*bg = */0.f);
    mEmitter->denseFill(CoordBBox(minEmtCoord, maxEmtCoord), /* value = */ 2.0, /*active = */ true);
    mEmitter->setTransform(mXform);
    mEmitter->setName("emitter");

    createDirichletVelocity2();

    mDensityCurr = FloatGrid::create(/*bg = */0.f);
    mDensityCurr->denseFill(CoordBBox(minBBoxIntrCoord, maxBBoxIntrCoord), /* value = */ 0.f, /* active = */ true);
    mDensityCurr->setTransform(mXform);
    mDensityCurr->setName("density_curr");
    mDensityCurr->topologyUnion(*mCollocatedMaskGrid);
    mDensityCurr->topologyIntersection(*mCollocatedMaskGrid);

    mVCurr = Vec3SGrid::create(/*bg = */Vec3s(0.f, 0.f, 0.f));
    mVCurr->denseFill(CoordBBox(minBBoxIntrCoord, mMaxStaggered), /* value = */ Vec3s(1.f, 0.f, 0.f), /* active = */ true);
    mVCurr->setTransform(mXform);
    mVCurr->setName("vel_curr");
    mVCurr->setGridClass(GRID_STAGGERED);
    mVCurr->tree().voxelizeActiveTiles();

    updateEmitter();
    applyDirichletVelocity(*mVCurr, -1);
}


void
SmokeSolver::initialize2() {
    std::cout << "initialize 2" << std::endl;
    using BBox = math::BBox<Vec3s>;

    mXform = math::Transform::createLinearTransform(mVoxelSize);
    float const padding = 1.f * mVoxelSize;
    mPadding = padding;
    float const centerY = 3.f;
    auto minBBox = Vec3s(0.f, 0.f, 0.f);
    auto maxBBox = Vec3s(5.f * mVoxelSize, 4.f * mVoxelSize, 4.f * mVoxelSize);
    // auto maxBBox = Vec3s(14.f, 6.f, 6.f);
    Coord minBBoxIntrCoord = mXform->worldToIndexNodeCentered(minBBox);
    Coord maxBBoxIntrCoord = mXform->worldToIndexNodeCentered(maxBBox);
    mMin = minBBoxIntrCoord;
    mMax = maxBBoxIntrCoord;
    mMaxStaggered = mMax + Coord(1);
    std::cout << "mMax = " << mMax << "\tmMaxStaggered = " << mMaxStaggered << std::endl;
    mCollocatedMaskGrid = BoolGrid::create(false /* background */);
    mCollocatedMaskGrid->denseFill(CoordBBox(minBBoxIntrCoord, maxBBoxIntrCoord), /* value = */ true, /* active = */ true);
    mCollocatedMaskGrid->setTransform(mXform);
    mCollocatedMaskGrid->setName("domain_mask");

    // Create an emitter and an emitter velocity
    auto minEmtW = Vec3s(0.f, 2.5f, 2.5f);
    auto maxEmtW = Vec3s(2 * mVoxelSize, 3.5f, 3.5f);
    Coord minEmtCoord = mXform->worldToIndexNodeCentered(minEmtW);
    Coord maxEmtCoord = mXform->worldToIndexNodeCentered(maxEmtW);
    mEmitter = FloatGrid::create(/*bg = */0.f);
    mEmitter->denseFill(CoordBBox(minEmtCoord, maxEmtCoord), /* value = */ 2.0, /*active = */ true);
    mEmitter->setTransform(mXform);
    mEmitter->setName("emitter");

    createDirichletVelocity2();

    mDensityCurr = FloatGrid::create(/*bg = */0.f);
    mDensityCurr->denseFill(CoordBBox(minBBoxIntrCoord, maxBBoxIntrCoord), /* value = */ 0.f, /* active = */ true);
    mDensityCurr->setTransform(mXform);
    mDensityCurr->setName("density_curr");
    mDensityCurr->topologyUnion(*mCollocatedMaskGrid);
    mDensityCurr->topologyIntersection(*mCollocatedMaskGrid);

    mVCurr = Vec3SGrid::create(/*bg = */Vec3s(0.f, 0.f, 0.f));
    mVCurr->denseFill(CoordBBox(minBBoxIntrCoord, mMaxStaggered), /* value = */ Vec3s(1.f, 0.f, 0.f), /* active = */ true);
    mVCurr->setTransform(mXform);
    mVCurr->setName("vel_curr");
    mVCurr->setGridClass(GRID_STAGGERED);
    mVCurr->tree().voxelizeActiveTiles();

    std::ostringstream ostr;
    ostr << "debug_initialize2.vdb";
    std::cerr << "\tWriting " << ostr.str() << std::endl;
    openvdb::io::File file(ostr.str());
    openvdb::GridPtrVec grids;
    grids.push_back(mVCurr);
    grids.push_back(mDensityCurr);
    file.write(grids);
    file.close();

    updateEmitter();
    applyDirichletVelocity(*mVCurr, -1);
}

void
SmokeSolver::addGravity(float const dt) {
    tree::LeafManager<Vec3STree> lm(mVCurr->tree());
    SmokeSolver::ApplyGravityOp op(dt, mGravity, mFlags);
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

    mDensityNext = advection.advect<FloatGrid, BoolGrid, SamplerT>(*mDensityCurr, *mCollocatedMaskGrid, dt);
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
    //         mDensityNext = a.advect<FloatGrid, BoolGrid, SamplerT>(*mDensityCurr, *mCollocatedMaskGrid, 0.1f);
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

    mVNext = advection.advect<Vec3SGrid, BoolGrid, SamplerT>(*mVCurr, *mCollocatedMaskGrid, dt);
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
    grids.push_back(mCollocatedMaskGrid);
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
    // 4 pm
    mDivBefore->topologyIntersection(*mDensityCurr);
    mDivBefore->setName("div_before");

    MaskGridType* domainMaskGrid = new MaskGridType(*mDivBefore); // match input grid's topology
    // 4 pm
    domainMaskGrid->topologyIntersection(*mDensityCurr);
    // tools::erodeActiveValues(domainMaskGrid->tree(), /*iterations=*/1, tools::NN_FACE, tools::IGNORE_TILES);
    // domainMaskGrid->topologyDifference(*mDirichletPressure);
    // domainMaskGrid->topologyDifference(*mDirichletVelocity);

    float divBefore = 0.f;
    auto divBeforeAcc = mDivBefore->getAccessor();
    for (auto iter = mDivBefore->beginValueOn(); iter; ++iter) {
        math::Coord ijk = iter.getCoord();
        auto val = divBeforeAcc.getValue(ijk);
        if (std::abs(val) > std::abs(divBefore)) {
            divBefore = val;
        }
    }
    std::cout << "\t== divergence before pp = " << divBefore << std::endl;

    math::pcg::State state = math::pcg::terminationDefaults<ValueType>();
    state.iterations = 100000;
    state.relativeError = state.absoluteError = epsilon;
    SmokeSolver::BoundaryOp bop(mDirichletVelocity, mDirichletPressure, mVoxelSize);
    util::NullInterrupter interrupter;
    FloatTree::Ptr fluidPressure = tools::poisson::solveWithBoundaryConditionsAndPreconditioner<PCT>(
        mDivBefore->tree(), domainMaskGrid->tree(), bop, state, interrupter, /*staggered=*/true);


    std::cout << "Projection Success: " << state.success << "\n";
    std::cout << "Iterations: " << state.iterations << "\n";
    std::cout << "Relative error: " << state.relativeError << "\n";
    std::cout << "Absolute error: " << state.absoluteError << "\n";

    // Note: need to dilate in order to do one-sided difference
    // because we use a staggered grid velocity field.
    FloatGrid::Ptr fluidPressureGrid = FloatGrid::create(fluidPressure);
    tools::dilateActiveValues(*fluidPressure, /*iterations=*/1, tools::NN_FACE, tools::IGNORE_TILES);

    fluidPressureGrid->setTransform(mXform);
    mPressure = fluidPressureGrid->copy();
    mPressure->setName("pressure");

    auto vCurrAcc = mVCurr->getAccessor();
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

    applyDirichletVelocity(*mVCurr, -2);
    mDivAfter = tools::divergence(*mVCurr);
    mDivAfter->setName("div_after");
    (mDivAfter->tree()).topologyIntersection(domainMaskGrid->tree());
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

    std::ostringstream ostr;
    ostr << "debug_divergence.vdb";
    std::cerr << "\tWriting " << ostr.str() << std::endl;
    openvdb::io::File file(ostr.str());
    openvdb::GridPtrVec grids;
    grids.push_back(mDivBefore);
    grids.push_back(mDivAfter);
    grids.push_back(mPressure);
    grids.push_back(mVCurr);
    grids.push_back(mDensityCurr);
    file.write(grids);
    file.close();

    //exit(0);



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
        updateEmitter();
        // addGravity(dt);
        pressureProjection4(true /* print */);
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
    grids.push_back(mCollocatedMaskGrid);

    file.write(grids);
    file.close();
}


void
SmokeSolver::writeVDBsDebug(int const frame) {
    std::cout << "Write VDBs Debug" << std::endl;
    std::ostringstream ss;
    ss << "INIT_DEBUG" << std::setw(3) << std::setfill('0') << frame << ".vdb";
    std::string fileName(ss.str());
    io::File file(fileName.c_str());

    openvdb::GridPtrVec grids;
    grids.push_back(mFlags);
    grids.push_back(mInteriorPressure);
    grids.push_back(mDirichletVelocity);
    grids.push_back(mVCurr);
    grids.push_back(mDensityCurr);
    grids.push_back(mDivBefore);
    grids.push_back(mPressure);

    file.write(grids);
    file.close();
}
}