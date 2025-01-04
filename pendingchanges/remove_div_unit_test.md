# Building
```shell
cmake .. -DOPENVDB_BUILD_UNITTESTS=ON --DOPENVDB_BUILD_EXAMPLES=ON
```

# Notes
```cpp
void
FlipSolver::substep(float const dt) {
    particlesToGrid();
    gridVelocityUpdate(dt);
    gridToParticles();
    updateParticles(dt);
}
```

`FlipSolver::gridVelocityUpdate(float const dt)`
- `velocityBCCorrection(mVCurr);`
   - this is enforcing the velocity boundary condition
- `pressureProjection(true /* print */);`
   - this is the guts of the algorithm
- `velocityBCCorrection(mVNext);`
   - also enforcing the velocity boundary condition
- `computeDivergence(mDivAfter, mVNext, "after");`

What are the grids that are needed?
- interior pressure. Look at its construction in `particlesToGrid`.

What's the problem you are trying to solve?
- You are given a velocity field and you want to make it divergence free subject to a certain boundary condition.