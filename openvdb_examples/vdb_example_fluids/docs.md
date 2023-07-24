# How to use the Poisson Solve in VDB

The important function to look at is: solveWithBoundaryConditionsAndPreconditioner.

The input of the algorithm is an input tree, a domain mask, a boundary operator, a state, an interrupter, and a boolean state to signify whether the velocity is staggered or not. Essentially, when the domain mask is not given, we use the same input tree as the domain mask.

Create an index tree from the domain mask. Create a mask of the interior voxels of the input tree by eroding the index tree.

If the domain mask is not given, then we use intree.

## Resources

OpenVDB:
https://www.openvdb.org/documentation/doxygen/codeExamples.html

OpenVDB point:
https://www.openvdb.org/documentation/doxygen/codeExamples.html#openvdbPointsHelloWorld

## Key

The key thing to understand are (1) the use of domain mask to tell the solver what dofs are considered fluid and (2) how to use a staggered grid velocity field in OpenVDB.


## How to create a grid

### Don't do this in production, do this in unit tests

## Working with MAC Grid
- Important to do:
  mVCurr->setGridClass(GRID_STAGGERED);
  before the divergence computation
- in the rasterization step, you need:
    TreeBase::Ptr tree =
        points::rasterizeTrilinear</*staggered=*/true, Vec3s>(points->tree(), "velocity");
## Pressure projection
- Need to erodeActiveVoxel to deal with Mac Grid during the divergence computation.
- Double check EXPAND_TILES or tools::IGNORE_TILES

## Topology:
 - Pressure interior dof is derived from doing erodeActiveVoxel on velocity grid.
 - Pressure full dof is derived either from mVCurr or erodeActiveVoxel followed by dilateActiveVoxel on the velocity grid.

 ## Where do I cut corners?
  - In the creation of fullPressure grid from fluidPressure

### Code and Answer

```
void
checkPoisson() {
    using TreeType = FloatTree;
    using ValueType = TreeType::ValueType;
    using std::sin;

    const ValueType zero = zeroVal<ValueType>();
    const double epsilon = math::Delta<ValueType>::value();

    int const N = 3;
    float const voxelSize = 1.0f/static_cast<float>(N);
    auto const xform = math::Transform::createLinearTransform(voxelSize);

    FloatTree::Ptr source(new FloatTree(0.f));
    source->fill(CoordBBox(Coord(1, 1, 1), Coord(N-1, N-1, N-1)), /* value = */0.f);
    FloatGrid::Ptr sourceGrid = Grid<FloatTree>::create(source);
    sourceGrid->setTransform(xform);
    auto srcAcc = sourceGrid->getAccessor();

    FloatTree::Ptr trueSln(new FloatTree(0.f));
    trueSln->fill(CoordBBox(Coord(0, 0, 0), Coord(N, N, N)), /* value = */0.f);
    FloatGrid::Ptr slnGrid = Grid<FloatTree>::create(trueSln);
    slnGrid->setTransform(xform);
    auto slnAcc = slnGrid->getAccessor();

    for (auto iter = sourceGrid->beginValueOn(); iter; ++iter) {
        auto ijk = iter.getCoord();
        auto xyz = xform->indexToWorld(ijk);
        // float sln = sin(2 * M_PI * xyz[0]) + sin(2 * M_PI * xyz[1]) + sin(2 * M_PI * xyz[2]);
        float sln = sin(2 * M_PI * xyz[0]);
        float rhs = -4.f * M_PI * M_PI * sln;
        srcAcc.setValue(ijk, rhs);
        slnAcc.setValue(ijk, sln);
        std::cout << "=== ijk" << ijk << " xyz = " << xyz << " = " << srcAcc.getValue(ijk) << std::endl;
    }

    SimpleExampleBoundaryOp bcOp(voxelSize);
    math::pcg::State state = math::pcg::terminationDefaults<ValueType>();
    state.iterations = 1000;
    state.relativeError = state.absoluteError = epsilon * 0.00000001;
    util::NullInterrupter interrupter;
    FloatTree::Ptr fluidPressure = tools::poisson::solveWithBoundaryConditions(
        sourceGrid->tree(), openvdb::tools::poisson::DirichletBoundaryOp<double>(), state, interrupter, /*staggered=*/true);

    std::cout << "Success: " << state.success << std::endl;
    std::cout << "Iterations: " << state.iterations << std::endl;
    std::cout << "Relative error: " << state.relativeError << std::endl;
    std::cout << "Absolute error: " << state.absoluteError << std::endl;
    std::cout << "before dilate solution->activeVoxelCount() =  " << fluidPressure->activeVoxelCount() << std::endl;
    std::cout << "epsilon = " << epsilon << std::endl;

    FloatGrid::Ptr fluidPressureGrid = FloatGrid::create(fluidPressure);
    auto numAcc = fluidPressureGrid->getAccessor();
    for (auto iter = sourceGrid->beginValueOn(); iter; ++iter) {
        auto ijk = iter.getCoord();
        auto xyz = xform->indexToWorld(ijk);
        float err = slnAcc.getValue(ijk) - numAcc.getValue(ijk) * voxelSize * voxelSize;
        float vsMult = numAcc.getValue(ijk) * voxelSize;
        float vsSqrMult = numAcc.getValue(ijk) * voxelSize * voxelSize;
        if (std::abs(err) > 1.0e-5) {
            std::cout << "ijk = " << ijk << " xyz = " << xyz << " true sln = " << slnAcc.getValue(ijk) << " ovdb sln = " << vsSqrMult << " err = " << err << std::endl;
        }
    }
}
```

Answer:
```
Success: 1
Iterations: 5
Relative error: 3.61712e-21
Absolute error: 1.23667e-19
before dilate solution->activeVoxelCount() =  8
epsilon = 1e-05
ijk = [1, 1, 1] xyz = [0.333333, 0.333333, 0.333333] true sln = 0.866025 ovdb sln = 0.759762 err = 0.106263
ijk = [1, 1, 2] xyz = [0.333333, 0.333333, 0.666667] true sln = 0.866025 ovdb sln = 0.759762 err = 0.106263
ijk = [1, 2, 1] xyz = [0.333333, 0.666667, 0.333333] true sln = 0.866025 ovdb sln = 0.759762 err = 0.106263
ijk = [1, 2, 2] xyz = [0.333333, 0.666667, 0.666667] true sln = 0.866025 ovdb sln = 0.759762 err = 0.106263
ijk = [2, 1, 1] xyz = [0.666667, 0.333333, 0.333333] true sln = -0.866025 ovdb sln = -0.759763 err = -0.106263
ijk = [2, 1, 2] xyz = [0.666667, 0.333333, 0.666667] true sln = -0.866025 ovdb sln = -0.759763 err = -0.106263
ijk = [2, 2, 1] xyz = [0.666667, 0.666667, 0.333333] true sln = -0.866025 ovdb sln = -0.759763 err = -0.106263
ijk = [2, 2, 2] xyz = [0.666667, 0.666667, 0.666667] true sln = -0.866025 ovdb sln = -0.759763 err = -0.106263
```