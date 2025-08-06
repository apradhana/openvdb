# Header Files:
There are 55 header files:
- Activate.h
- ChangeBackground.h
- Clip.h
- Composite.h
- Count.h
- Dense.h
- DenseSparseTools.h
- Diagnostics.h
- FastSweeping.h
- Filter.h
- FindActiveValues.h
- GridOperators.h
- GridTransformer.h
- Interpolation.h
- LevelSetAdvect.h
- LevelSetDilatedMesh.h
- LevelSetFilter.h
- LevelSetFracture.h
- LevelSetMeasure.h
- LevelSetMorph.h
- LevelSetPlatonic.h
- LevelSetRebuild.h
- LevelSetSphere.h
- LevelSetTracker.h
- LevelSetTubes.h
- LevelSetUtil.h
- Mask.h
- Merge.h
- MeshToVolume.h
- Morphology.h
- MultiResGrid.h
- NodeVisitor.h
- ParticleAtlas.h
- ParticlesToLevelSet.h
- PointAdvect.h
- PointIndexGrid.h
- PointPartitioner.h
- PointScatter.h
- PointsToMask.h
- PoissonSolver.h
- PotentialFlow.h
- Prune.h
- RayIntersector.h
- RayTracer.h
- SignedFloodFill.h
- Statistics.h
- tools.md
- tools.txt
- TopologyToLevelSet.h
- ValueTransformer.h
- VectorTransformer.h
- VelocityFields.h
- VolumeAdvect.h
- VolumeToMesh.h
- VolumeToSpheres.h

# More function

- [ ] Activate.h [2]: changing the topology of a grid.
  - [ ] activate
        Purpose:
        Marks as active any inactive tiles or voxels in the given grid or tree whose values are equal to a specified value (optionally within a given tolerance).
        Effect: change the topology of a grid.
        Signature:
        ```cpp
        template<typename GridOrTree>
        void activate(
            GridOrTree&,
            const typename GridOrTree::ValueType& value,
            const typename GridOrTree::ValueType& tolerance = zeroVal<typename GridOrTree::ValueType>(),
            const bool threaded = true
        );
        ```
        Example:
        ```cpp
        openvdb::FloatGrid::Ptr grid = openvdb::FloatGrid::create(0.0f);
        // Activate all voxels with value exactly 0.0f
        openvdb::tools::activate(*grid, 0.0f);
        ```
  - [ ] deactivate
        Purpose:
        Marks as inactive any active tiles or voxels in the given grid or tree whose values are equal to a specified value (optionally within a given tolerance).
        Effect: change the topology of a grid.
        Signature:
        ```cpp
        template<typename GridOrTree>
        void deactivate(
            GridOrTree&,
            const typename GridOrTree::ValueType& value,
            const typename GridOrTree::ValueType& tolerance = zeroVal<typename GridOrTree::ValueType>(),
            const bool threaded = true
        );
        ```
        Example:
        ```cpp
        openvdb::FloatGrid::Ptr grid = openvdb::FloatGrid::create(1.0f);
        // Deactivate all voxels with value exactly 1.0f
        openvdb::tools::deactivate(*grid, 1.0f);
        ```

- [ ] ChangeBackground.h
- [ ] Clip.h
- [ ] Composite.h
- [ ] Count.h
- [ ] Dense.h
- [ ] DenseSparseTools.h
- [ ] Diagnostics.h
- [ ] FastSweeping.h
- [ ] Filter.h
- [ ] FindActiveValues.h
- [ ] GridOperators.h
- [ ] GridTransformer.h
- [ ] Interpolation.h
- [ ] LevelSetAdvect.h
- [ ] LevelSetDilatedMesh.h
- [ ] LevelSetFilter.h
- [ ] LevelSetFracture.h
- [ ] LevelSetMeasure.h
- [ ] LevelSetMorph.h
- [ ] LevelSetPlatonic.h
- [ ] LevelSetRebuild.h
- [ ] LevelSetSphere.h
- [ ] LevelSetTracker.h
- [ ] LevelSetTubes.h
- [ ] LevelSetUtil.h
- [ ] Mask.h
- [ ] Merge.h
- [ ] MeshToVolume.h
- [ ] Morphology.h
- [ ] MultiResGrid.h
- [ ] NodeVisitor.h
- [ ] ParticleAtlas.h
- [ ] ParticlesToLevelSet.h
- [ ] PointAdvect.h
- [ ] PointIndexGrid.h
- [ ] PointPartitioner.h
- [ ] PointScatter.h
- [ ] PointsToMask.h
- [ ] PoissonSolver.h
- [ ] PotentialFlow.h
- [ ] Prune.h
- [ ] RayIntersector.h
- [ ] RayTracer.h
- [ ] SignedFloodFill.h
- [ ] Statistics.h
- [ ] tools.txt
- [ ] TopologyToLevelSet.h
- [ ] ValueTransformer.h
- [ ] VectorTransformer.h
- [ ] VelocityFields.h
- [ ] VolumeAdvect.h
- [ ] VolumeToMesh.h
- [ ] VolumeToSpheres.h

# Prompt

Give me a list of the public function API in this file and tell me what is the purpose of each. Give a tiny example on how to use the function. Give it in a list form marked with - [ ] in markdown format.

