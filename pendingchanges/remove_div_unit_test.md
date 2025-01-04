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
- A little bit of complication comes from the fact that the boundary condition is formulated in terms of pressure interior cells, while the velocity are the faces of these.

What is a Neumann pressure cell?

What is a Dirichlet pressure cell?

If the B.C. is all Neumann, we will have null-space.

Perhaps the most crucial part of the implementation is how a user defines a "Boundary Operator" functor that gets passed into the class. The user defined boundary operator is called for each dof in the mask where its neighbor lies outside of the mask. This boundary operator is what defines whether a dof is a Dirichlet or a Neumann pressure, and depending on the choice, the user needs to tell the solver how to modify the diagonal entry of the Laplacian matrix and the right hand side corresponding to that dof. 

You can have a flag grid. This is to be fed into the BoundaryOperator, which is then fed into tools::poisson::solveWithBoundaryConditionsAndPreconditioner.

