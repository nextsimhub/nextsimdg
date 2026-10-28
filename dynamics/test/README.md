# Testcases for the Dynamics Core

## Model Array Test Cases

There are three tests checking that the ModelArray framework works correctly:

- CGModelArray_test tests the conversion between the dynamics code CGVector and ModelArray.
- DGModelArray_test tests the conversion between the dynamics code DGVector and ModelArray.
- ParametricMesh_test tests the processing of ModelArray in the ParametricMesh class, as well as the reading of `.smesh` files.

## Linear Advection Test Cases

These test cases test the DG transport scheme for dG0, dG1 and dG2

### Advection_test

A smooth initial function is transported for times [0,2pi] in the domain [0,1]^2 using the rotational velocity field v = (y-1/2, x-1/2)

We print out the final mass (for testing mass conservation) and the error after one revolution (numerical solution compared to initial solution).


### AdvectionPeriodicBC_test

Here we test the implementation of the boundary handling: A smooth initial density f(x,y) = exp(-50 (x-1/2)^2 - 50 (y-1/2)^2 ) is transported to the upper right corner then to the lower left and back to the middle.

When 'hitting' the corners, parts of the smooth function are transported out of the domain. This tests 'free outflow behavior'. When transporting back, the boundary value is zero.

We check the mass at final time:

- We compare it to the exact value (integrating f(x,y)) over the remaining part for measuring convergence
- We compare it to previously stored reference values for automatic checking of failure.

The three test cases are using different domains / discretizations:

- Uniform discretization of [0,1] x [0,1] consisting of N x N elements with mesh size h x h
- Uniform discretization of [0,2] x [0,1] consisting of 2N x N elements with mesh size h x h
- Non-uniform discretization of [0,2] x [0,1] consisting of N x N elements using mesh size (2h) x h
