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

What is a Neumann pressure cell?  What is a Dirichlet pressure cell?  
// 0 Neumann
// 1 interior
// 4 dirichlet pressure. In this setup it's on the right. It means that it's not a collider on the right.
// Neumann pressure means Dirichlet velocity.
If the B.C. is all Neumann, we will have null-space.

The most crucial part of the implementation is how a user defines a "Boundary Operator" functor that gets passed into the class. The user defined boundary operator is called for each dof in the mask where its neighbor lies outside of the mask. This boundary operator is what defines whether a dof is a Dirichlet or a Neumann pressure, and depending on the choice, the user needs to tell the solver how to modify the diagonal entry of the Laplacian matrix and the right hand side corresponding to that dof. 

You can have a flag grid. This is to be fed into the BoundaryOperator, which is then fed into tools::poisson::solveWithBoundaryConditionsAndPreconditioner. More specifically, the BoundaryOperator also needs dirichletVelocity grid to enforce the boundary condition correctly.


# Debugging the pressure projection in the smoke solver

## Regression Test
```shell
Writing flags.vdb
create dirichlet velocity 4
Write VDBs Debug
frame = 0 substep = 10
update emitter
done with update emitter
apply dirichlet velocity begins
pressure projection 4
== divergence before pp = -0.0833333
Projection Success: 0
Iterations: 82
Relative error: 0.0297418
Absolute error: 0.894733
apply dirichlet velocity begins
== divergence after pp = -0.894734
```
## QnAs

**Q:** Do I need to call applyDirichletVelocity again?

## Steps
[X] Pass regression test
[X] Add computeDivergence and computeLInfinity APIs
[ ] Add v(x, y, z) = (x^2, y^2, z^2)

```shell
Writing flags.vdb
create dirichlet velocity 4
Write VDBs Debug
frame = 0 substep = 10
update emitter
done with update emitter
apply dirichlet velocity begins
pressure projection debug
Divergence before = 30.0833
== divergence before pp = 30.0833
Projection Success: 0
Iterations: 82
Relative error: 0.0297418
Absolute error: 0.894733
apply dirichlet velocity begins
== divergence after pp = -0.894734
```

```shell
Writing flags.vdb
create dirichlet velocity 4
Write VDBs Debug
frame = 0 substep = 10
update emitter
done with update emitter
apply dirichlet velocity begins
pressure projection debug
apply dirichlet velocity begins
Divergence before = -74.6
Projection Success: 0
Iterations: 49
Relative error: 1.41231
Absolute error: 105.358
apply dirichlet velocity begins
Divergence after = 105.358
```

With FLIP:
```shell
==frame = 2==
Divergence before = -1.44118
Projection success: 1
Projection iterations: 162
Projection relative error: 7.63079e-06
Projection absolute error: 1.09973e-05
Pressure->activeVoxelCount() =  45097
Divergence after = 1.09427e-05
```

# Branch
```
https://github.com/apradhana/openvdb/tree/test
```