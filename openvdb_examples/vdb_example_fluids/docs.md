# How to use the Poisson Solve in VDB

The important function to look at is: solveWithBoundaryConditionsAndPreconditioner.

The input of the algorithm is an input tree, a domain mask, a boundary operator, a state, an interrupter, and a boolean state to signify whether the velocity is staggered or not. Essentially, when the domain mask is not given, we use the same input tree as the domain mask.

Create an index tree from the domain mask. Create a mask of the interior voxels of the input tree by eroding the index tree.

If the domain mask is not given, then we use intree.

## Key

The key thing to understand are (1) the use of domain mask to tell the solver what dofs are considered fluid and (2) how to use a staggered grid velocity field in OpenVDB.

## How to create a grid

### Don't do this in production, do this in unit tests
