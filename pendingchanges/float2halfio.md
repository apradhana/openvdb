[ ] Create half grid
[ ] Read from stream2memory

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