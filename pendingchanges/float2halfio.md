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

## Dec 5
Call stack from readGrid to Archive::doReadGrid:
1. TestGridIO::testConvertFloatToHalf (TestGridIO.cc:362)
   └─ baseGridHalf = fileHalf.readGrid(nameIter.gridName());

2. File::readGrid(const Name&) (File.cc:602-605)
   └─ return readGridByName(name, BBoxd());

3. File::readGridByName(const Name&, const BBoxd&) (File.cc:617)
   └─ grid = (clip ? readGrid(gd, bbox) : readGrid(gd));
      [assuming clip=false, calls readGrid(gd)]

4. File::readGrid(const GridDescriptor&) const (File.cc:827-830)
   └─ return Impl::readGrid(*this, gd, Impl::NoBBox());

5. File::Impl::readGrid<NoBBox>(const File&, const GridDescriptor&, const NoBBox&) (File.cc:47-58)
   └─ GridBase::Ptr grid = file.createGrid(gd);
   └─ gd.seekToGrid(file.inputStream());
   └─ unarchive(file, grid, gd, bbox);

6. File::Impl::unarchive(const File&, GridBase::Ptr&, const GridDescriptor&, NoBBox) (File.cc:61-66)
   └─ file.Archive::readGrid(grid, gd, file.inputStream());

7. Archive::readGrid(GridBase::Ptr, const GridDescriptor&, std::istream&, ScalarConversion) (Archive.cc:1271-1280)
   └─ readGridCompression(is);
   └─ doReadGrid(grid, gd, is, NoBBox());

8. doReadGrid<NoBBox>(GridBase::Ptr, const GridDescriptor&, std::istream&, const NoBBox&) (Archive.cc:1137-1265)
   [Final destination - performs the actual grid reading]
   
This tells me that I only need to add ScalarConversion argument inside File::Impl::unarchive

## Dec 2
Archive::readGrid()
  → RootNode::readTopologyWithValueType<SourceValueT>()
    → InternalNode::readTopologyWithValueType<SourceValueT>()
      → LeafNode::readTopologyWithValueType<SourceValueT>()

Grid::readBuffers(is, CoordBBox) [Grid.h:1652]
  → Tree::readBuffers(is, CoordBBox, saveFloatAsHalf) [Tree.h:1313]
    → RootNode::readBuffers(is, CoordBBox, fromHalf) [RootNode.h:2522]
      → InternalNode::readBuffers(is, CoordBBox, fromHalf) [InternalNode.h:3264]
        → LeafNode::readBuffers(is, CoordBBox, fromHalf) [LeafNode.h:1398]

From PointIndexLeafNode::readBuffers(CoordBBox) (line 1739 in tools/PointIndexGrid.h)
PointIndexLeafNode<T, Log2Dim>::readBuffers(std::istream& is, const CoordBBox& bbox, bool fromHalf)
{
    // Read and clip voxel values.
    BaseLeaf::readBuffers(is, bbox, fromHalf);
    ...
}

Seems like you need to add readBuffersWithValue.

## Dec 1

Why we can't mix template with virtual function.

Option 2: Use Type Erasure with a Base Reader Class
If you need polymorphic behavior, create a non-template base reader interface:

// In io/Compression.h or appropriate header
namespace io {
    class ReaderBase {
    public:
        virtual ~ReaderBase() = default;
        // Add virtual methods needed for reading
    };

    template<typename ValueT, typename ValueFromT>
    class TypedConvertingReader : public ReaderBase {
        // Your implementation
    };
}

// Then in GridBase:
virtual void readTopology(std::istream&, io::ReaderBase& reader) = 0;


Explanation of the Error
The error "templates may not be 'virtual'" is a fundamental C++ language restriction. You cannot combine the virtual and template keywords on the same function. Here's why:
The Core Problem
Virtual Functions (Runtime Polymorphism):
Use a vtable (virtual function table) created at compile time
Function calls are resolved at runtime through pointer/reference indirection
The vtable has fixed entries - the compiler needs to know exactly which functions exist
Template Functions (Compile-time Polymorphism):
Are instantiated on-demand when you use them with specific types
The compiler generates a new function for each unique set of template arguments
This happens at compile time, potentially creating many different functions
Why They're Incompatible