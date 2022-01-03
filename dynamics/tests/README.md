# Testcases for the Dynamics Core

## Linear Advection Test Cases

These test cases test the dg transport scheme for dG0, dG1 and dG2

### example1

A smooth initial function is transported for times [0,2pi] in the domain [0,1]^2 using the rotational velocity field v = (y-1/2, x-1/2)

We print out the final mass (for testing mass conservation) and the error after one revolution (numerical solution compared to initial solution). 

FAILURE CHECKING STILL MISSING

### example2a, example2b, example 2c

Here we test the implementation of the boundary handling: A smooth initial density f(x,y) = exp(-50 (x-1/2)^2 - 50 (y-1/2)^2 ) is transported to the upper right corner then to the lower left and back to the middle. 

When 'hitting' the corners, parts of the smooth function are transported out of the domain. This tests 'free outflow behavior'. When transporting back, the boundary value is zero. 

We check the mass at final time:

- We compare it to the exact value (integrating f(x,y)) over the remaining part for measuring convergence
- We compare it to previously stored reference values for automatic checking of failure. 

The three test cases are using different domains / discretizations:

- example2a uniform discretization of [0,1] x [0,1] consisting of N x N elements with mesh size h x h
- example2b uniform discretization of [0,2] x [0,1] consisting of 2N x N elements with mesh size h x h
- example2c non-uniform discretization of [0,2] x [0,1] consisting of N x N elements using mesh size (2h) x h

### multithreading\_square / multithreading\_rectangle1/2

As example1 testing the transport in a rotational velocity field. We check the convergence of the scheme for increased mesh fineness. 

- multithreading\_square uniform discretization of [0,1] x [0,1] consisting of N x N elements with mesh size h x h
- multithreading\_rectangle1 uniform discretization of [0,2] x [0,1] consisting of 2N x N elements with mesh size h x h
- multithreading\_rectangle2 non-uniform discretization of  [0,2] x [0,1] consisting of N x N elements using mesh size (2h) x h

FAILURE CHECKING NOT IMPLEMENTED SO FAR


## Test cases for the momentum equation

