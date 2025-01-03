# Building
```shell
cmake .. -DOPENVDB_BUILD_UNITTESTS=ON --DOPENVDB_BUILD_EXAMPLES=ON
```

# Notes
`FlipSolver::gridVelocityUpdate(float const dt)`
- `velocityBCCorrection(mVCurr);`
   - this is enforcing the velocity boundary condition
- `pressureProjection(true /* print */);`
   - this is the guts of the algorithm
- `velocityBCCorrection(mVNext);`
   - also enforcing the velocity boundary condition
- `computeDivergence(mDivAfter, mVNext, "after");`