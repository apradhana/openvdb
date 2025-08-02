// Copyright Contributors to the OpenVDB Project
// SPDX-License-Identifier: Apache-2.0
//
/// @file main.cc
///
/// @brief Simple ray tracer for OpenVDB volumes
///
/// @note This is intended mainly as an example of how to ray-trace
/// OpenVDB volumes.  It is not a production-quality renderer.

#include <openvdb/openvdb.h>
#include <openvdb/openvdb.h>
#include <openvdb/math/Mat.h>
#include <openvdb/tools/LevelSetSphere.h>
#include <openvdb/tools/LevelSetPlatonic.h>
#include <openvdb/tools/MeshToVolume.h>
#include <fstream>
#include <iostream>
#include <vector>


namespace internal {
struct CopyValuesOp
{
    CopyValuesOp(openvdb::HalfGrid::Ptr hg,
                 openvdb::FloatGrid::Ptr fg) :
                 hg(hg),
                 fg(fg) {}

    template <typename T>
    void operator()(T &node, size_t) const
    {
        auto halfAcc = hg->getAccessor();
        for (typename T::ValueAllIter iter = node.beginValueAll(); iter; ++iter) {
            auto ijk = iter.getCoord();
            iter.setValue(halfAcc.getValue(ijk));
        }
    }

    openvdb::HalfGrid::Ptr hg;
    openvdb::FloatGrid::Ptr fg;
};// CopyValuesOp
} // namespace internal

template<typename GridT>
struct TestResult {
    double duration;
    typename GridT::Ptr grid;
    int activeVoxelCount;

    TestResult(double duration, typename GridT::Ptr grid, int activeVoxelCount) : duration(duration), grid(grid), activeVoxelCount(activeVoxelCount) {}

    void print() const {
        std::cout << "duration = " << duration << " seconds" << std::endl;
    }
};

template<typename HalfGridT>
void
convertHalfToFloatGrid(typename HalfGridT::Ptr hg, const openvdb::FloatGrid::Ptr fg, std::string gridName)
{
    using namespace openvdb;

    fg->setName(gridName);
    fg->setTransform(hg->transform().copy());
    fg->tree().topologyUnion(hg->tree());

    tree::LeafManager<FloatTree> lm(fg->tree());
    ::internal::CopyValuesOp op(hg, fg);
    lm.foreach(op);

    auto fAcc = fg->getAccessor();
    auto hAcc = hg->getAccessor();
    float maxDif = 0.f;
    for (auto iter = fg->beginValueOn(); iter; ++iter) {
        math::Coord const ijk = iter.getCoord();
        auto const fv = fAcc.getValue(ijk);
        auto const hv = hAcc.getValue(ijk);
        float const dif = std::abs(fv - hv);
        if (dif > maxDif) {
            maxDif = dif;
        }
    }

    if (hg->getGridClass() == GRID_LEVEL_SET) {
        fg->setGridClass(GRID_LEVEL_SET);
        openvdb::tools::changeLevelSetBackground(fg->tree(), 3.f);
    }
    std::cout << "convertHalfToFloatGrid::maxDif = " << maxDif << "\tactiveVoxelCount dif = " << (int)(fg->activeVoxelCount() - hg->activeVoxelCount()) << "\n";
}

void
saveTestResults(TestResult<openvdb::FloatGrid> floatResult, TestResult<openvdb::HalfGrid> halfResult, std::string fileName)
{
    using namespace openvdb;

    FloatGrid::Ptr halfInFloat = FloatGrid::create();
    convertHalfToFloatGrid<HalfGrid>(gridPtrCast<HalfGrid>(halfResult.grid), halfInFloat, halfResult.grid->getName() + "_in_float");

    GridPtrVec grids;
    grids.push_back(floatResult.grid);
    grids.push_back(halfInFloat);

    openvdb::io::File file(fileName);
    file.write(grids);
    file.close();
}


template<typename GridType>
TestResult<GridType>
testLevelSetSphereImpl(float radius, float voxelSize, float halfWidth, openvdb::Vec3f center, std::string name)
{
    const tbb::tick_count start = tbb::tick_count::now();
    typename GridType::Ptr grid = openvdb::tools::createLevelSetSphere<GridType>(
        radius, center, voxelSize, halfWidth);
    grid->setName(name);
    grid->setGridClass(openvdb::GRID_LEVEL_SET);
    const tbb::tick_count end = tbb::tick_count::now();
    const double duration = (end - start).seconds();

    return TestResult<GridType>(duration, grid, grid->activeVoxelCount());
}

void testLevelSetSphere() {
    const std::string fileName = "testLevelSetSphere.vdb";
    const std::string logTitle = "==== Test create level set sphere ====\n";

    std::cout << logTitle;

    const float radius = 50.0f;
    const openvdb::Vec3f center(0.0f, 0.0f, 0.0f);
    const float voxelSize = 0.1f;
    const float halfWidth = 3.0f; // narrow band half-width in voxels

    TestResult<openvdb::HalfGrid> halfResult = testLevelSetSphereImpl<openvdb::HalfGrid>(radius, voxelSize, halfWidth, center, "half_sphere");
    TestResult<openvdb::FloatGrid> floatResult = testLevelSetSphereImpl<openvdb::FloatGrid>(radius, voxelSize, halfWidth, center, "float_sphere");

    openvdb::GridPtrVec grids;
    grids.push_back(halfResult.grid);
    grids.push_back(floatResult.grid);
    openvdb::io::File file(fileName);
    file.write(grids);
    file.close();

    halfResult.print();
    floatResult.print();
}


template<typename GridType>
TestResult<GridType>
testLevelSetPlatonicImpl(int faceCount, float scale, openvdb::Vec3f center, float voxelSize, float halfWidth, std::string name)
{
    const tbb::tick_count start = tbb::tick_count::now();
    typename GridType::Ptr grid = openvdb::tools::createLevelSetPlatonic<GridType>(
        faceCount, scale, center, voxelSize, halfWidth);
    grid->setName(name);
    grid->setGridClass(openvdb::GRID_LEVEL_SET);
    const tbb::tick_count end = tbb::tick_count::now();
    const double duration = (end - start).seconds();

    return TestResult<GridType>(duration, grid, grid->activeVoxelCount());
}

void testLevelSetPlatonic()
{
    const std::string fileName = "testLevelSetPlatonic.vdb";
    const std::string logTitle = "==== Test create level set platonic ====\n";

    std::cout << logTitle;

    const int faceCount = 8; // 4=Tetrahedron, 6=Cube, 8=Octahedron, 12=Dodecahedron, 20=Icosahedron
    const float scale = 30.0f;
    const openvdb::Vec3f center(0.0f, 0.0f, 0.0f);
    const float voxelSize = 0.1f;
    const float halfWidth = 3.0f; // narrow band half-width in voxels

    // Create a level set platonic solid as a FloatGrid
    auto floatResult = testLevelSetPlatonicImpl<openvdb::FloatGrid>(faceCount, scale, center, voxelSize, halfWidth, "float_octahedron");
    auto halfResult = testLevelSetPlatonicImpl<openvdb::HalfGrid>(faceCount, scale, center, voxelSize, halfWidth, "half_octahedron");

    openvdb::GridPtrVec grids;
    grids.push_back(floatResult.grid);
    openvdb::io::File file(fileName);
    file.write(grids);
    file.close();

    floatResult.print();
    halfResult.print();
}


template<typename GridType>
TestResult<GridType>
testMeshToVolumeImpl(std::vector<openvdb::Vec3s> points, std::vector<openvdb::Vec3I> triangles, std::vector<openvdb::Vec4I> quads, float exteriorBandWidth, float interiorBandWidth, openvdb::math::Transform::Ptr transform, std::string name)
{
    using namespace openvdb;

    typename GridType::Ptr grid = nullptr;

    const tbb::tick_count start = tbb::tick_count::now();
    if (!triangles.empty()) {
        grid = tools::meshToLevelSet<GridType>(*transform, points, triangles, exteriorBandWidth);
    } else if (!quads.empty()) {
        grid = tools::meshToLevelSet<GridType>(*transform, points, quads, exteriorBandWidth);
    }
    const tbb::tick_count end = tbb::tick_count::now();
    const double duration = (end - start).seconds();
    if (!grid) {
        std::cerr << "Failed to create VDB grid from mesh." << std::endl;
    }

    return TestResult<GridType>(duration, grid, grid->activeVoxelCount());
}

void parseObjFile(std::string objFile, std::vector<openvdb::Vec3s>& points, std::vector<openvdb::Vec3I>& triangles, std::vector<openvdb::Vec4I>& quads)
{
    using namespace openvdb;

    std::ifstream in(objFile);
    if (!in) {
        std::cerr << "Failed to open OBJ file: " << objFile << std::endl;
    }
    const tbb::tick_count start = tbb::tick_count::now();
    std::string line;
    while (std::getline(in, line)) {
        std::istringstream iss(line);
        std::string type;
        iss >> type;
        if (type == "v") {
            float x, y, z;
            iss >> x >> y >> z;
            points.emplace_back(x, y, z);
        } else if (type == "f") {
            std::vector<int> idx;
            std::string vert;
            while (iss >> vert) {
                std::istringstream viss(vert);
                int i;
                viss >> i;
                idx.push_back(i - 1); // OBJ is 1-based
            }
            if (idx.size() == 3) {
                triangles.emplace_back(idx[0], idx[1], idx[2]);
            } else if (idx.size() == 4) {
                quads.emplace_back(idx[0], idx[1], idx[2], idx[3]);
            }
        }
    }

    if (points.empty() || (triangles.empty() && quads.empty())) {
        std::cerr << "OBJ file does not contain valid mesh data." << std::endl;
    }
    const tbb::tick_count end = tbb::tick_count::now();
    const double duration = (end - start).seconds();
    std::cout << "OBJ file parsing duration = " << duration << " seconds" << std::endl;
}


void testMeshToVolume()
{
    using namespace openvdb;

    std::string objFile = "/home/andre/Desktop/dragon.obj";
    std::string vdbFile = "/home/andre/Desktop/dragon_meshtovolume.vdb" ;
    std::string logTitle = "==== Test mesh to volume ====\n";

    std::cout << logTitle;

    // Parse OBJ
    std::vector<Vec3s> points;
    std::vector<Vec3I> triangles;
    std::vector<Vec4I> quads;
    parseObjFile(objFile, points, triangles, quads);

    // Create a transform (identity, voxel size = 1.0)
    float voxelSize = 0.2f;
    math::Transform::Ptr transform = math::Transform::createLinearTransform(voxelSize);

    // Convert mesh to level set
    float exteriorBandWidth = 3.0f, interiorBandWidth = 3.0f;
    auto floatResult = testMeshToVolumeImpl<openvdb::FloatGrid>(points, triangles, quads, exteriorBandWidth, interiorBandWidth, transform, "float_dragon");
    auto halfResult = testMeshToVolumeImpl<openvdb::HalfGrid>(points, triangles, quads, exteriorBandWidth, interiorBandWidth, transform, "half_dragon");

    saveTestResults(floatResult, halfResult, vdbFile);

    floatResult.print();
    halfResult.print();
}

int main()
{
    openvdb::initialize();

    testLevelSetSphere();
    testLevelSetPlatonic();
    testMeshToVolume();

    return 0;
}
