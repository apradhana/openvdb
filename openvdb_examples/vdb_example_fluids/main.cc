// Copyright Contributors to the OpenVDB Project
// SPDX-License-Identifier: MPL-2.0

#include <algorithm>
#include <fstream>
#include <iostream>
#include <math.h>
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
#include <openvdb/tools/GridOperators.h> // for divergence
using namespace openvdb;

class Vector3 {
public:
    float x, y, z;

    Vector3(float _x, float _y, float _z) : x(_x), y(_y), z(_z) {}
};

class SmokeSolver {
public:
    SmokeSolver();

    void initialize();
    void foobar();

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
    std::cout << "foobar begins" << std::endl;
    HeatMethod hm(0.1f);
    hm.createFlags();
    hm.createDirichletBC();
    std::cout << "foobar ends" << std::endl;
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
    convertToOnesAndZeros();
    SmokeSolver solver;
    solver.foobar();
    // solver.render();
    // testPoissonSolve();
    // testDivergence();
}
