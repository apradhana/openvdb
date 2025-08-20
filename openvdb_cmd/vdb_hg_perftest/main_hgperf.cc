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
#include <openvdb/tools/VolumeToMesh.h>
#include <openvdb/tools/ChangeBackground.h>
#include <openvdb/tools/LevelSetMeasure.h>
#include <openvdb/tools/LevelSetAdvect.h>
#include <openvdb/tools/LevelSetTracker.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <sstream>
#include <iomanip>

namespace internal {
template<typename InputGridT, typename OutputGridT>
struct ConvertValuesOp
{
    ConvertValuesOp(typename InputGridT::Ptr inputGrid,
                    typename OutputGridT::Ptr outputGrid) :
                    inputGrid(inputGrid),
                    outputGrid(outputGrid) {}

    template <typename T>
    void operator()(T &node, size_t) const
    {
        auto inputAcc = inputGrid->getAccessor();
        for (typename T::ValueAllIter iter = node.beginValueAll(); iter; ++iter) {
            auto ijk = iter.getCoord();
            iter.setValue(inputAcc.getValue(ijk));
        }
    }

    typename InputGridT::Ptr inputGrid;
    typename OutputGridT::Ptr outputGrid;
};// ConvertValuesOp

template<typename VecGridT>
struct RotationVelocityFieldOp
{
    using ValueT = typename VecGridT::ValueType;

    RotationVelocityFieldOp(typename VecGridT::Ptr vel) : velGrid(vel) {}

    template <typename T>
    void operator()(T &node, size_t) const
    {
        auto velAcc = velGrid->getAccessor();
        auto xform = velGrid->transform();
        for (typename T::ValueAllIter iter = node.beginValueAll(); iter; ++iter) {
            auto ijk = iter.getCoord();
            auto xyz = xform.indexToWorld(ijk);
            ValueT newVel(0.0f, 0.0f, 0.0f);
            iter.setValue(newVel);
        }
    }

    typename VecGridT::Ptr velGrid;
};// RotationVelocityFieldOp
} // namespace internal

template<typename GridT>
struct TestResult {
    double duration;
    typename GridT::Ptr grid;
    int activeVoxelCount;

    TestResult(double duration, typename GridT::Ptr grid, int activeVoxelCount) : duration(duration), grid(grid), activeVoxelCount(activeVoxelCount) {}

    void print(std::string prefix) const {
        std::cout << prefix << " duration = " << duration << " seconds, activeVoxelCount = " << activeVoxelCount << "\n";
    }
};

template<typename InputGridT, typename OutputGridT>
void
convertPayload(typename InputGridT::Ptr inputGrid, typename OutputGridT::Ptr outputGrid, std::string gridName)
{
    using namespace openvdb;

    using OutputTreeType = typename OutputGridT::TreeType;
    using ValueType = typename OutputGridT::ValueType;

    outputGrid->setName(gridName);
    outputGrid->setTransform(inputGrid->transform().copy());
    outputGrid->tree().topologyUnion(inputGrid->tree());

    tree::LeafManager<OutputTreeType> lm(outputGrid->tree());
    ::internal::ConvertValuesOp<InputGridT, OutputGridT> op(inputGrid, outputGrid);
    lm.foreach(op);

    auto outAcc = outputGrid->getAccessor();
    auto inAcc = inputGrid->getAccessor();
    float maxDif = 0.f;
    for (auto iter = outputGrid->beginValueOn(); iter; ++iter) {
        math::Coord const ijk = iter.getCoord();
        auto const outv = outAcc.getValue(ijk);
        auto const inv = inAcc.getValue(ijk);
        float const dif = std::abs(outv - inv);
        if (dif > maxDif) {
            maxDif = dif;
        }
    }
    std::cout << "convertPayload::maxDif = " << maxDif << "\tactiveVoxelCount dif = " << (int)(outputGrid->activeVoxelCount() - inputGrid->activeVoxelCount()) << "\n";


    if (inputGrid->getGridClass() == GRID_LEVEL_SET) {
        outputGrid->setGridClass(GRID_LEVEL_SET);
        openvdb::tools::changeLevelSetBackground(outputGrid->tree(), ValueType(3.f));
    }
}

template<typename HalfGridT>
void
convertHalfToFloatGrid(typename HalfGridT::Ptr hg, const openvdb::FloatGrid::Ptr fg, std::string gridName)
{
    using namespace openvdb;

    fg->setName(gridName);
    fg->setTransform(hg->transform().copy());
    fg->tree().topologyUnion(hg->tree());

    tree::LeafManager<FloatTree> lm(fg->tree());
    ::internal::ConvertValuesOp<HalfGrid, FloatGrid> op(hg, fg);
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
    convertPayload<HalfGrid, FloatGrid>(gridPtrCast<HalfGrid>(halfResult.grid), halfInFloat, halfResult.grid->getName() + "_in_float");

    GridPtrVec grids;
    grids.push_back(floatResult.grid);
    grids.push_back(halfInFloat);

    io::File file(fileName);
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

    saveTestResults(floatResult, halfResult, fileName);

    halfResult.print("half_sphere ");
    floatResult.print("float_sphere");
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

    saveTestResults(floatResult, halfResult, fileName);

    floatResult.print("float_octahedron");
    halfResult.print("half_octahedron ");
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

    floatResult.print("float_dragon");
    halfResult.print("half_dragon ");
}

void _writeObjFile(const std::string& filename,
                  const std::vector<openvdb::Vec3s>& points,
                  const std::vector<openvdb::Vec3I>& triangles,
                  const std::vector<openvdb::Vec4I>& quads)
{
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open " << filename << " for writing" << std::endl;
        return;
    }

    // Write vertices
    for (const auto& point : points) {
        file << "v " << point[0] << " " << point[1] << " " << point[2] << std::endl;
    }

    // Write triangular faces
    for (const auto& triangle : triangles) {
        file << "f " << (triangle[0] + 1) << " " << (triangle[1] + 1) << " " << (triangle[2] + 1) << std::endl;
    }

    // Write quad faces
    for (const auto& quad : quads) {
        file << "f " << (quad[0] + 1) << " " << (quad[1] + 1) << " " << (quad[2] + 1) << " " << (quad[3] + 1) << std::endl;
    }

    file.close();
    std::cout << "Mesh written to " << filename << std::endl;
}

void testVolumeToMesh()
{
    using namespace openvdb;
    std::string vdbFile = "dragon.vdb";
    std::string objFileFloat = "dragon_mesh_float.obj";
    std::string objFileHalf = "dragon_mesh_half.obj";
    std::string logTitle = "==== Test volumeToMesh ====";
    std::cout << logTitle << std::endl;

    // Open the VDB file
    io::File file(vdbFile);
    file.open();
    GridPtrVecPtr grids = file.getGrids();
    if (!grids || grids->empty()) {
        std::cerr << "No grids found in " << vdbFile << std::endl;
        return;
    }
    // Find the first FloatGrid
    FloatGrid::Ptr floatGridDragon = nullptr;
    for (auto& baseGrid : *grids) {
        floatGridDragon = gridPtrCast<FloatGrid>(baseGrid);
        if (floatGridDragon) break;
    }
    if (!floatGridDragon) {
        std::cerr << "No FloatGrid found in " << vdbFile << std::endl;
        return;
    }
    HalfGrid::Ptr halfGridDragon = HalfGrid::create();
    convertPayload<FloatGrid, HalfGrid>(floatGridDragon, halfGridDragon, floatGridDragon->getName() + "_half");
    std::cout << "floatGridDragon activeVoxelCount: " << floatGridDragon->activeVoxelCount() << " halfGridDragon activeVoxelCount: " << halfGridDragon->activeVoxelCount() << std::endl;

    // Run volumeToMesh
    std::vector<openvdb::Vec3s> pointsFloat, pointsHalf;
    std::vector<openvdb::Vec3I> trianglesFloat, trianglesHalf;
    std::vector<openvdb::Vec4I> quadsFloat, quadsHalf;
    double isovalue = 0.0;
    double adaptivity = 0.0;
    bool relaxDisorientedTriangles = true;
    double durationFloat, durationHalf;
    {
        const tbb::tick_count start = tbb::tick_count::now();
        openvdb::tools::volumeToMesh(*floatGridDragon, pointsFloat, trianglesFloat, quadsFloat, isovalue, adaptivity, relaxDisorientedTriangles);
        const tbb::tick_count end = tbb::tick_count::now();
        durationFloat = (end - start).seconds();
    }
    {
        const tbb::tick_count start = tbb::tick_count::now();
        openvdb::tools::volumeToMesh(*halfGridDragon, pointsHalf, trianglesHalf, quadsHalf, isovalue, adaptivity, relaxDisorientedTriangles);
        const tbb::tick_count end = tbb::tick_count::now();
        durationHalf = (end - start).seconds();
    }

    _writeObjFile(objFileFloat, pointsFloat, trianglesFloat, quadsFloat);
    _writeObjFile(objFileHalf, pointsHalf, trianglesHalf, quadsHalf);

    std::cout << "volumeToMesh results for float grid:" << std::endl;
    std::cout << "  duration:  " << durationFloat << " seconds" << std::endl;
    std::cout << "  mesh saved to:  " << objFileFloat << std::endl;
    std::cout << "  points:    " << pointsFloat.size() << std::endl;
    std::cout << "  triangles: " << trianglesFloat.size() << std::endl;
    std::cout << "  quads:     " << quadsFloat.size() << std::endl;
    std::cout << "volumeToMesh results for half grid:" << std::endl;
    std::cout << "  duration:  " << durationHalf << " seconds" << std::endl;
    std::cout << "  mesh saved to:  " << objFileHalf << std::endl;
    std::cout << "  points:    " << pointsHalf.size() << std::endl;
    std::cout << "  triangles: " << trianglesHalf.size() << std::endl;
    std::cout << "  quads:     " << quadsHalf.size() << std::endl;
}

void testDirectConversionToHalfGrid() {
    using namespace openvdb;
    std::string fileName = "./dragon.vdb";

    io::File file(fileName);
    file.open(false /* delay load*/, io::MappedFile::Notifier(), io::Archive::ScalarConversion::Half);

    HalfGrid::Ptr grid;
    // Loop over the names of all of the grids in the file.
    for (io::File::NameIterator nameIter = file.beginName();
        nameIter != file.endName(); ++nameIter)
    {
        std::string gridName = nameIter.gridName();
        grid = gridPtrCast<HalfGrid>(file.readGrid(gridName));
        std::cout << "gridName = " << gridName << "grid = " << grid << std::endl;
    }
}

void testChangeLevelSetBackground()
{
    using namespace openvdb;
    std::string vdbFile = "dragon.vdb";
    std::string logTitle = "==== Test changeLevelSetBackground ====";
    std::cout << logTitle << std::endl;

    io::File file(vdbFile);
    file.open();
    GridPtrVecPtr grids = file.getGrids();
    if (!grids || grids->empty()) {
        std::cerr << "No grids found in " << vdbFile << std::endl;
        return;
    }
    FloatGrid::Ptr grid = nullptr;
    for (auto& baseGrid : *grids) {
        grid = gridPtrCast<FloatGrid>(baseGrid);
        if (grid) break;
    }
    if (!grid) {
        std::cerr << "No FloatGrid found in " << vdbFile << std::endl;
        return;
    }

    HalfGrid::Ptr halfGrid = HalfGrid::create();
    auto halfGridBgBefore = halfGrid->background();
    halfGrid->setName(grid->getName() + "_half");
    halfGrid->setTransform(grid->transform().copy());
    halfGrid->tree().topologyUnion(grid->tree());

    openvdb::tools::changeLevelSetBackground(halfGrid->tree(), 3.f);

    std::cout << "changeLevelSetBackground results:" << std::endl;
    std::cout << "  background before = " << halfGridBgBefore << ", after = " << halfGrid->background() << std::endl;
}

void testLevelSetMeasure()
{
    using namespace openvdb;
    std::string vdbFile = "dragon.vdb";
    std::string logTitle = "==== Test LevelSetMeasure ====";
    std::cout << logTitle << std::endl;

    // Open the VDB file
    io::File file(vdbFile);
    file.open();
    GridPtrVecPtr grids = file.getGrids();
    if (!grids || grids->empty()) {
        std::cerr << "No grids found in " << vdbFile << std::endl;
        return;
    }

    // Find the first FloatGrid
    FloatGrid::Ptr floatGrid = nullptr;
    for (auto& baseGrid : *grids) {
        floatGrid = gridPtrCast<FloatGrid>(baseGrid);
        if (floatGrid) break;
    }
    if (!floatGrid) {
        std::cerr << "No FloatGrid found in " << vdbFile << std::endl;
        return;
    }

    // Create HalfGrid version
    HalfGrid::Ptr halfGrid = HalfGrid::create();
    convertPayload<FloatGrid, HalfGrid>(floatGrid, halfGrid, floatGrid->getName() + "_half");

    std::cout << "Float grid name: " << floatGrid->getName() << " half grid name: " << halfGrid->getName() << std::endl;
    std::cout << "FloatGrid activeVoxelCount: " << floatGrid->activeVoxelCount() << " half grid activeVoxelCount: " << halfGrid->activeVoxelCount() << std::endl;

    // Test LevelSetMeasure functions for FloatGrid
    try {
        Real floatArea = tools::levelSetArea(*floatGrid, true);
        Real floatVolume = tools::levelSetVolume(*floatGrid, true);
        int floatEuler = tools::levelSetEulerCharacteristic(*floatGrid);
        int floatGenus = tools::levelSetGenus(*floatGrid);
        Real halfArea = tools::levelSetArea(*halfGrid, true);
        Real halfVolume = tools::levelSetVolume(*halfGrid, true);
        int halfEuler = tools::levelSetEulerCharacteristic(*halfGrid);
        int halfGenus = tools::levelSetGenus(*halfGrid);

        std::cout << "Measure area, volume, Euler characteristic, genus separately:" << std::endl;
        std::cout << "  Diff Surface Area: " << floatArea - halfArea << " world units squared" << std::endl;
        std::cout << "  Diff Volume: " << floatVolume - halfVolume << " world units cubed" << std::endl;
        std::cout << "  Diff Euler Characteristic: " << floatEuler - halfEuler << std::endl;
        std::cout << "  Diff Genus: " << floatGenus - halfGenus << std::endl << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "FloatGrid and HalfGrid measurement error: " << e.what() << std::endl;
    }

    // Test LevelSetMeasure class directly for more detailed analysis
    try {
        tools::LevelSetMeasure<FloatGrid> measureFloat(*floatGrid);
        tools::LevelSetMeasure<HalfGrid> measureHalf(*halfGrid);

        Real areaFloat = measureFloat.area(true);
        Real volumeFloat = measureFloat.volume(true);
        Real avgMeanCurvatureFloat = measureFloat.avgMeanCurvature(true);
        Real avgGaussianCurvatureFloat = measureFloat.avgGaussianCurvature(true);
        int eulerCharFloat = measureFloat.eulerCharacteristic();
        int genusFloat = measureFloat.genus();

        Real areaHalf = measureHalf.area(true);
        Real volumeHalf = measureHalf.volume(true);
        Real avgMeanCurvatureHalf = measureHalf.avgMeanCurvature(true);
        Real avgGaussianCurvatureHalf = measureHalf.avgGaussianCurvature(true);
        int eulerCharHalf = measureHalf.eulerCharacteristic();
        int genusHalf = measureHalf.genus();

        std::cout << "Measure area, volume, Euler characteristic, genus using LevelSetMeasure class" << std::endl;
        std::cout << "  Diff Surface Area: " << areaFloat - areaHalf << " world units squared" << std::endl;
        std::cout << "  Diff Volume: " << volumeFloat - volumeHalf << " world units cubed" << std::endl;
        std::cout << "  Diff Average Mean Curvature: " << avgMeanCurvatureFloat - avgMeanCurvatureHalf << std::endl;
        std::cout << "  Diff Average Gaussian Curvature: " << avgGaussianCurvatureFloat - avgGaussianCurvatureHalf << std::endl;
        std::cout << "  Diff Euler Characteristic: " << eulerCharFloat - eulerCharHalf << std::endl;
        std::cout << "  Diff Genus: " << genusFloat - genusHalf << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "LevelSetMeasure class measurement error: " << e.what() << std::endl;
    }
}
void testLevelSetAdvection()
{
    using namespace openvdb;
    std::string logTitle = "==== Test LevelSetAdvection ====";
    std::cout << logTitle << std::endl;

    const int faceCount = 4; // 4=Tetrahedron, 6=Cube, 8=Octahedron, 12=Dodecahedron, 20=Icosahedron
    const float scale = 1.0f;
    const openvdb::Vec3f center(0.0f, 0.0f, 0.0f);
    const float voxelSize = 0.1f;
    const float halfWidth = 3.0f; // narrow band half-width in voxels

    // Create a level set platonic solid as a FloatGrid
    auto floatLS = testLevelSetPlatonicImpl<openvdb::FloatGrid>(faceCount, scale, center, voxelSize, halfWidth, "float_tetrahedron");
    auto halfLS = testLevelSetPlatonicImpl<openvdb::HalfGrid>(faceCount, scale, center, voxelSize, halfWidth, "half_tetrahedron");

    // Create velocity field
    auto cubeRes = testLevelSetPlatonicImpl<openvdb::FloatGrid>(6, scale, center, voxelSize, halfWidth, "float_cube");
    auto cubeLS = cubeRes.grid;
    openvdb::Vec3fGrid::Ptr velocityField = openvdb::Vec3fGrid::create(openvdb::Vec3f(0.0f, 0.0f, 0.0f));
    velocityField->setTransform(cubeLS->transform().copy());
    velocityField->tree().topologyUnion(cubeLS->tree());
    velocityField->setName("velocity_field");

    tree::LeafManager<Vec3fGrid::TreeType> lm(velocityField->tree());
    ::internal::RotationVelocityFieldOp<Vec3fGrid> op(velocityField);
    lm.foreach(op);

}

int main()
{
    openvdb::initialize();

    // authoring level set
    testLevelSetSphere();
    testLevelSetPlatonic();

    // level set operations
    testChangeLevelSetBackground();

    // level set measurement
    testLevelSetMeasure();

    // conversion from mesh to level set
    testMeshToVolume();
    testVolumeToMesh();

    testDirectConversionToHalfGrid();

    return 0;
}
