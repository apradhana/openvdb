// Copyright Contributors to the OpenVDB Project
// SPDX-License-Identifier: MPL-2.0

#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <openvdb/openvdb.h>
#include <openvdb/points/PointConversion.h>
#include <openvdb/points/PointCount.h>
#include <openvdb/util/logging.h>
#include <openvdb/tree/NodeManager.h> // for post processing bool grid 

class Vector3 {
public:
    float x, y, z;

    Vector3(float _x, float _y, float _z) : x(_x), y(_y), z(_z) {}
};

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
// struct TrueNearZeroKernel {
//     TrueNearZeroKernel(float epsilon) : mEpsilon(epsilon) {}
// 
//     // Root node
//     void operator()(typename openvdb::FloatTree::RootNodeType& node, size_t = 1) const {
//         for (auto iter = node.beginValueAll(); iter; ++iter) {
//             if (std::abs(*iter) < mEpsilon) {
//                 iter.
//             }
//         }
//     }
// 
//     // Internal nodes
//     template<typename NodeT>
//     void operator()(NodeT& node, size_t = 1) const
//     {
//         for (auto iter = node.beginValueAll(); iter; ++iter) {
//             if (*iter == -std::numeric_limits<SdfValueT>::max()) {
//                 iter.setValue(mMin);
//             }
//             if (*iter == std::numeric_limits<SdfValueT>::max()) {
//                 iter.setValue(mMax);
//             }
//         }
//     }
// 
//     // Leaf nodes
//     void operator()(typename SdfTreeT::LeafNodeType& leaf, size_t = 1) const
//     {
//         for (auto iter = leaf.beginValueOn(); iter; ++iter) {
//             if (*iter == -std::numeric_limits<SdfValueT>::max()) {
//                 iter.setValue(mMin);
//             }
//             if (*iter == std::numeric_limits<SdfValueT>::max()) {
//                 iter.setValue(mMax);
//             }
//         }
//     }
//   float mEpsilon;
// };// FastSweeping::TrueNearZeroKernel


void convertToBool() {
    openvdb::io::File file("/home/andre/dev/openvdb_aswf/_build_fluids/openvdb_examples/vdb_example_fluids/c_skin_mid.vdb");
    file.open();
    openvdb::GridBase::Ptr baseGrid;
    openvdb::io::File::NameIterator nameIter = file.beginName();
    baseGrid = file.readGrid(nameIter.gridName());
    file.close();
    openvdb::FloatGrid::Ptr sdf = openvdb::gridPtrCast<openvdb::FloatGrid>(baseGrid);

    openvdb::BoolTree bTree(sdf->tree(), false, openvdb::TopologyCopy());

    std::cout << "sdf->activeVoxelCount() = " << sdf->activeVoxelCount() << std::endl;
    auto xform = sdf->transform();
    xform.print();
}

void foobar() {
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

// TO BUILD:
// mkdir build
// cd build
// cmake -DOPENVDB_BUILD_EXAMPLES=ON -DOPENVDB_BUILD_VDB_EXAMPLE_FLUIDS=ON ../
// make -j 8

int
main(int argc, char *argv[])
{
    openvdb::initialize();

    foobar();
}
