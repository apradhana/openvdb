[X] Create half grid
[X] Read from stream2memory
[ ] Check whether adding ConverterReader in GridDescriptor changes the ABI.
[ ] Add an argument in readTopology in RootNode for ConvertingReader.
[ ] Create convertingReader in GridDescriptor.cc while we are creating grid, then pass it down to Archive::doRead and later to RootNode::readTopology and LeafNode::read

How to run test:
```shell
cd _build
```
## What mismatch do you have between the istream and the grid?
doReadGrid (Archive.cc:1137)
├─ readTopology (line 1215)
│  ├─ 1. RootNode background ............ +2 bytes FIRST ERROR
│  ├─ 2. RootNode tile values ........... +2 bytes per tile
│  └─ For each InternalNode:
│     ├─ 3a. Inactive value 0 ........... +2 bytes
│     ├─ 3b. Inactive value 1 ........... +2 bytes
│     └─ 3c. Tile buffer data ........... +2×count bytes
│
└─ readBuffers (line 1216)
   └─ For each LeafNode:
      ├─ 4a. Inactive value 0 ........... +2 bytes
      ├─ 4b. Inactive value 1 ........... +2 bytes
      └─ 4c. Voxel buffer data ........... +2×SIZE bytes (≈1024 bytes per leaf!)

# What test are you running with ./vdb_test?
```shell
./vdb_test --gtest_filter=TestLeafIOTest*:TestFile*:TestGridDescriptor*:TestTree*
```

What is a filtering step in the work by Ivo?
```shell
./vdb_test --gtest_filter=TestLeafIOTest*:TestFile*:TestGridDescriptor*:TestTree*:-TestLeafIOTest.testBufferInt:TestLeafIOTest.testBufferFloat:TestLeafIOTest.testBufferDouble:TestLeafIOTest.testBufferByte:TestLeafIOTest.testBufferVec3R:TestFile.testWriteGrid:TestFile.testWriteMultipleGrids:TestFile.testReadGridDescriptors:TestFile.testEmptyGridIO:TestFile.testDelayedLoadMetadata:TestGridDescriptor.testIO:TestTree.testHalf:TestTree.testIO
```

## What statistics are you checking?
RootNode::readTopologyWithValueType - sizeof(SourceValueT): 4 bytes
RootNode::readTopologyWithValueType - ValueType: float
RootNode::readTopologyWithValueType - sizeof(ValueType): 4 bytes
RootNode::readTopologyWithValueType - mBackground: 0.150024
RootNode::readTopologyWithValueType - numTiles: 0
RootNode::readTopologyWithValueType - numChildren: 8

## Other directory:
Original branch:
feature/hg_vdb_render_ghurst/feature/half_grid_support

## What cmake command did you use?
cmake
cmake .. -DCMAKE_BUILD_TYPE=Release -DOPENVDB_BUILD_UNITTESTS=ON

## The right answer for dragon level set:
=== HalfGrid Statistics ===
Background: 0.300049
Leaf nodes: 124166
Node counts by level: L0=124166 L1=318 L2=8 L3=1
Non-leaf nodes: 327
Root tiles: 0
Root child nodes: 8