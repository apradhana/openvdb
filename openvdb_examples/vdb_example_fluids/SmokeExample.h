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

    void initialize4();

    void substep(float const dt, int const frame);

    // Make the velocity on the grid to be divergence free
    void pressureProjection(bool print);

    void updateEmitter();
    void createDirichletVelocity();
    void createVCurr(bool print);
    void createDensityCurr();
    void createFlags(bool print);
    void createInteriorPressure();
    void applyDirichletVelocity4(Vec3SGrid& vecGrid, int frame);

    void swapGrids();

    void addGravity(float const dt);

    void advectDensity(float const dt);

    void advectVelocity(float const dt, int const frame);

    void writeVDBs(int const frame);

    void writeVDBsDebug(int const frame);

    struct BoundaryOp4 {
        BoundaryOp4(Int32Grid::ConstPtr flags,
                    Vec3SGrid::ConstPtr dirichletVelocity,
                    float const voxelSize) :
                    flags(flags),
                    dirichletVelocity(dirichletVelocity),
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
            } else if (isDirichletPressure) {
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
    Vec3s mGravity = Vec3s(0.f, -2.f, 0.f);
    Vec3s mPushVelocity = Vec3s(0.2f, 0.f, 0.f);
    math::Transform::Ptr mXform;
    float mPadding;
    Coord mMin;
    Coord mMax;
    Coord mMaxStaggered;

    Int32Grid::Ptr mFlags;
    FloatGrid::Ptr mDensityCurr;
    FloatGrid::Ptr mDensityNext;
    Vec3SGrid::Ptr mVCurr;
    Vec3SGrid::Ptr mVNext;

    BoolGrid::Ptr mInteriorPressure;

    FloatGrid::Ptr mSphere;
    FloatGrid::Ptr mEmitter;
    Vec3SGrid::Ptr mDirichletVelocity;

    FloatGrid::Ptr mDivBefore;
    FloatGrid::Ptr mDivAfter;
    FloatGrid::Ptr mPressure;
};


SmokeSolver::SmokeSolver(float const voxelSize) : mVoxelSize(voxelSize)
{
    initialize4();
}

 void
 SmokeSolver::createFlags(bool print)
 {
    mFlags = Int32Grid::create(/* bg = */ 0); // Neumann pressure
    mFlags->denseFill(CoordBBox(mMin, mMax), /* value = */ 1, /* active = */ true);
    mFlags->setTransform(mXform);
    mFlags->setName("flags");
    auto flagsAcc = mFlags->getAccessor();
    auto sphereAcc = mSphere->getAccessor();
    for (auto iter = mFlags->beginValueOn(); iter; ++iter) {
        math::Coord ijk = iter.getCoord();

        if (ijk[0] == mMax[0]) {
            flagsAcc.setValue(ijk, 4); // Dirichlet
        }
        if (ijk[0] == mMin[0] /* left face */ ||
            ijk[1] == mMin[1] /* bottom face */ ||
            ijk[1] == mMax[1] /* top face */ ||
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
SmokeSolver::createInteriorPressure()
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
 SmokeSolver::createDensityCurr()
 {
    mDensityCurr = FloatGrid::create(/*bg = */0.f);
    mDensityCurr->setTransform(mXform);
    mDensityCurr->setName("density_curr");
    mDensityCurr->denseFill(CoordBBox(mMin, mMax), /* value = */ 0.f, /* active = */ true);
    mDensityCurr->topologyUnion(*mFlags);
    mDensityCurr->topologyIntersection(*mInteriorPressure);
 }


 void
 SmokeSolver::createVCurr(bool print)
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

    auto flagsAcc = mFlags->getConstAccessor();
    auto velAcc = mVCurr->getAccessor();
    for (auto iter = mVCurr->beginValueOn(); iter; ++iter) {
        auto ijk = iter.getCoord();
        auto im1jk = ijk.offsetBy(-1, 0, 0);
        auto ijm1k = ijk.offsetBy(0, -1, 0);
        auto ijkm1 = ijk.offsetBy(0, 0, -1);

        if (flagsAcc.getValue(ijk) == 0) {
            if (flagsAcc.getValue(im1jk) == 0 &&
                flagsAcc.getValue(ijm1k) == 0 &&
                flagsAcc.getValue(ijkm1) == 0) {
                    velAcc.setValueOff(ijk);
                }
        }
    }

 }


 void
 SmokeSolver::initialize4()
 {
    std::cout << "initialize 4" << std::endl;
    using BBox = math::BBox<Vec3s>;
    mPushVelocity = Vec3s(3.f, 0.f, 0.f);
    mXform = math::Transform::createLinearTransform(mVoxelSize);
    float const padding = 1.f * mVoxelSize;
    mPadding = padding;
    float const centerY = 3.f;
    auto minBBox = Vec3s(0.f, 0.f, 0.f);
    // works:
    // auto maxBBox = Vec3s(3 * mVoxelSize, 3 * mVoxelSize, 3 * mVoxelSize);
    //auto maxBBox = Vec3s(1.f + mVoxelSize, 0.4f + mVoxelSize, 0.4f + mVoxelSize);
    auto maxBBox = Vec3s(5.f + mVoxelSize, 1.5f + mVoxelSize, 1.5f + mVoxelSize);
    Coord minBBoxIntrCoord = mXform->worldToIndexNodeCentered(minBBox);
    Coord maxBBoxIntrCoord = mXform->worldToIndexNodeCentered(maxBBox);
    mMin = minBBoxIntrCoord;
    mMax = maxBBoxIntrCoord;
    mMaxStaggered = mMax + Coord(1);

    const float radius = 0.2 * maxBBox[1];
    const openvdb::Vec3f center(2.5f / 7.f * maxBBox[0], 0.5f * maxBBox[1], 0.5f * maxBBox[2]);
    mSphere = tools::createLevelSetSphere<openvdb::FloatGrid>(radius, center, mVoxelSize, 2 /* width */);

    // Create an emitter and an emitter velocity
    auto minEmtW = Vec3s(0.f, 0.5f*maxBBox[1]-0.1*maxBBox[1], 0.5f*maxBBox[2]-0.1*maxBBox[2]);
    auto maxEmtW = Vec3s(10 * mVoxelSize, 0.5f*maxBBox[1]+0.2*maxBBox[1], 0.5f*maxBBox[2]+0.2*maxBBox[2]);
    Coord minEmtCoord = mXform->worldToIndexNodeCentered(minEmtW);
    Coord maxEmtCoord = mXform->worldToIndexNodeCentered(maxEmtW);
    mEmitter = FloatGrid::create(/*bg = */0.f);
    mEmitter->denseFill(CoordBBox(minEmtCoord, maxEmtCoord), /* value = */ 2.0, /*active = */ true);
    mEmitter->setTransform(mXform);
    mEmitter->setName("emitter");

    createFlags(false);
    createInteriorPressure();
    createVCurr(false);
    createDensityCurr();
    createDirichletVelocity();
    writeVDBsDebug(0);
 }

 void
 SmokeSolver::pressureProjection(bool print)
 {
    std::cout << "pressure projection 4" << std::endl;
    using TreeType = FloatTree;
    using ValueType = TreeType::ValueType;
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
    SmokeSolver::BoundaryOp4 bop(mFlags, mVCurr, mVoxelSize);
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


    // for (auto iter = mPressure->beginValueOn(); iter; ++iter) {
    //     auto ijk = iter.getCoord();
    //     std::cout << "p" << ijk << " = " << pressureAcc.getValue(ijk) << std::endl;
    // }



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

    // if (!state.success) {
    //     std::ostringstream ostr;
    //     ostr << "debug_velocity_fail.vdb";
    //     std::cerr << "\tWriting " << ostr.str() << std::endl;
    //     openvdb::io::File file(ostr.str());
    //     openvdb::GridPtrVec grids;
    //     grids.push_back(mVCurr);
    //     grids.push_back(mPressure);
    //     file.write(grids);
    //     exit(0);
    // }
 }

 void
 SmokeSolver::createDirichletVelocity() {
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
        if (ijk[0] == 1) {
            drcAcc.setValue(ijk, mPushVelocity);
        }
        // auto im1jk = ijk.offsetBy(-1, 0, 0);
        // auto ijm1k = ijk.offsetBy(0, -1, 0);
        // auto ijkm1 = ijk.offsetBy(0, 0, -1);

        // Vec3s val;

        // if (flagsAcc.getValue(ijk) == 0) {
        //     // is a full Neumann
        //     val = mPushVelocity;
        // } else {
        //     if(flagsAcc.getValue(im1jk) == 0) {
        //         // neighboring a Neumann pressure in the x face
        //         val[0] = mPushVelocity[0];
        //     }
        //     if(flagsAcc.getValue(ijm1k) == 0) {
        //         // neighboring a Neumann pressure in the y face
        //         val[1] = mPushVelocity[1];
        //     }
        //     if(flagsAcc.getValue(ijkm1) == 0) {
        //         // neighboring a Neumann pressure in the z face
        //         val[2] = mPushVelocity[2];
        //     }
        // }
        // drcAcc.setValue(ijk, val);
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
SmokeSolver::swapGrids()
{
    mDensityCurr = mDensityNext;
    mDensityCurr->setName("density_curr");
    mVCurr = mVNext;
    mVCurr->setName("vel_curr");
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
SmokeSolver::addGravity(float const dt) {
    tree::LeafManager<Vec3STree> lm(mVCurr->tree());
    SmokeSolver::ApplyGravityOp op(dt, mGravity, mFlags);
    lm.foreach(op);
}


void
SmokeSolver::advectDensity(float const dt)
{
    std::cout << "Advect Density" << std::endl;
    using AdvT = openvdb::tools::VolumeAdvection<Vec3fGrid, true /* Staggered */>;
    using SamplerT = openvdb::tools::Sampler<1>;

    AdvT advection(*mVCurr);
    advection.setIntegrator(tools::Scheme::SEMI);
    advection.setLimiter(tools::Scheme::REVERT);
    advection.setSubSteps(1);

    mDensityNext = advection.advect<FloatGrid, BoolGrid, SamplerT>(*mDensityCurr, *mInteriorPressure, dt);
    mDensityNext->setName("density_next");
}

void
SmokeSolver::advectVelocity(float const dt, const int frame)
{
    std::cout << "Advect Velocity" << std::endl;
    using AdvT = openvdb::tools::VolumeAdvection<Vec3SGrid, true /* staggered */>;
    using SamplerT = openvdb::tools::Sampler<1, true /* staggered */>;

    AdvT advection(*mVCurr);
    advection.setIntegrator(tools::Scheme::SEMI);
    advection.setSubSteps(1);

    mVNext = advection.advect<Vec3SGrid, BoolGrid, SamplerT>(*mVCurr, *mInteriorPressure, dt);
    // mVNext = advection.advect<Vec3SGrid, SamplerT>(*mVCurr, dt);
    mVNext->setGridClass(GRID_STAGGERED);
    mVNext->setName("vel_next");
}


void
SmokeSolver::substep(float const dt, int const frame) {
    updateEmitter();
    addGravity(dt);
    applyDirichletVelocity4(*mVCurr, -1);
    pressureProjection(false);
    advectVelocity(dt, frame);
    advectDensity(dt);
    swapGrids();
}


void
SmokeSolver::render() {
    float const dt = 1.f/24.f;
    for (int frame = 0; frame < 600; ++frame) {
        float const numSubStep = 50.f;
        for (int i = 0; i < static_cast<int>(numSubStep); ++i) {
             substep(dt / numSubStep, frame);
        }
        writeVDBsDebug(frame);
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
    grids.push_back(mDensityCurr);
    grids.push_back(mVCurr);
    grids.push_back(mPressure);

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
    grids.push_back(mEmitter);
    grids.push_back(mDirichletVelocity);
    grids.push_back(mFlags);
    grids.push_back(mInteriorPressure);
    grids.push_back(mDirichletVelocity);
    grids.push_back(mVCurr);
    grids.push_back(mDensityCurr);
    grids.push_back(mDivBefore);
    grids.push_back(mDivAfter);
    grids.push_back(mPressure);

    file.write(grids);
    file.close();
}
}