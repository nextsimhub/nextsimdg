"""
Defines the DG and CG basis functions used in the different code generation functions.

The functions are specific to the unit square [0,1]^2.
"""

import numpy as np


def dgdofs(d):
    """
    Compute the number of local unknowns in the DG space.

    :param d: Number of unknowns per element depending on gauss degree
    :return: the number of local unknowns in the dg spaces
    """
    if d==0:
        return 1
    elif d==1:
        return 3
    elif d==2:
        return 6
    else:
        msg = "dG3 and higher is not implemented"
        raise AssertionError(msg)

def cgdofs(d):
    """
    Compute the number of local unknowns in the CG space.

    :param d: Number of unknowns per element depending on gauss degree
    :return: the number of local unknowns in the cg spaces
    """
    if d==1:
        return 4
    elif d==2:
        return 9
    else:
        msg = "only cg1 and cg2 are supported"
        raise AssertionError(msg)

def dgbasis(j,x,y):
    """
    Evaluate the dG-basis functions on [0,1]^2 in (x,y).

    dG 0:  1
    dG 1:  x, y
    dG 2:  x^2, y^2, xy
    dG X:  x^2y, xy^2.
    """
    if j==0:
        return 1.
    elif j==1:
        return x-0.5
    elif j==2:
        return y-0.5
    elif j==3:
        return (x-0.5)*(x-0.5)-1.0/12.0
    elif j==4:
        return (y-0.5)*(y-0.5)-1.0/12.0
    elif j==5:
        return (x-0.5)*(y-0.5)
    elif j==6:
        return (y-0.5)*((x-0.5)*(x-0.5)-1.0/12.0)
    elif j==7:
        return (x-0.5)*((y-0.5)*(y-0.5)-1.0/12.0)
    else:
        print("dG3 and higher not implemented (yet)")
        raise AssertionError

def dx_dgbasis(j,x,y):
    if j==0:
        return 0.
    elif j==1:
        return 1.
    elif j==2:
        return 0.
    elif j==3:
        return 2.0*(x-0.5)
    elif j==4:
        return 0.
    elif j==5:
        return (y-0.5)
    elif j==6:
        return (y-0.5)*(2.*(x-0.5))
    elif j==7:
        return ((y-0.5)*(y-0.5)-1.0/12.0)
    else:
        print("dG3 and higher not implemented (yet)")
        raise AssertionError

def dy_dgbasis(j,x,y):
    if j in (0, 1):
        return 0.
    elif j==2:
        return 1.
    elif j==3:
        return 0.
    elif j==4:
        return 2.*(y-0.5)
    elif j==5:
        return (x-0.5)
    elif j==6:
        return ((x-0.5)*(x-0.5)-1.0/12.0)
    elif j==7:
        return (x-0.5)*(2.*(y-0.5))
    else:
        print("dG3 and higher not implemented (yet)")
        raise AssertionError



def dgbasis_edge(j,x):
    """Evaluate 1d dG-basis on the edge [0,1]."""
    if j==0:
        return 1.
    elif j==1:
        return x-0.5
    elif j==2:
        return (x-0.5)*(x-0.5)-1.0/12.0
    else:
        print("dG3 and higher not implemented (yet)")
        raise AssertionError

# Inverse element mass matrix for the dg methods
inversemass = np.array([1., 12., 12., 180., 180., 144., 2160., 2160.])


def CGbasis1d(cg,j,x):
    """Evaluate the CG(2)-basis functions on [0,1]^2 in (x,y)."""
    if cg==1:
        if j==0:
            return 1.0-x
        elif j==1:
            return x
        else:
            print("CG1basis1d only for j=0,1")
            raise AssertionError

    elif cg==2:
        if j==0:
            return 2.0*(x-0.5)*(x-1.0)
        elif j==1:
            return 4.0*x*(1.0-x)
        elif j==2:
            return 2.0*x*(x-0.5)
        else:
            print("CGbasis1d only for j=0,1,2")
            raise AssertionError
    else:
        print("only CG1 CG2")
        raise AssertionError

def CGbasis1d_dX(cg,j,x):
    """... and its derivative."""
    if cg==1:
        if j==0:
            return -1
        else:
            return 1
    elif j==0:
        return 4.0*x-3.0
    elif j==1:
        return 4.0-8.0*x
    elif j==2:
        return 4.0*x-1.0
    else:
        print("CGbasis1d only for j=0,1,2")
        raise AssertionError

def CGbasisfunction(cg,j,x,y):
    jx = j%(cg+1)
    jy = j//(cg+1)
    return CGbasis1d(cg,jx,x)*CGbasis1d(cg,jy,y)

def CGbasisfunction_dX(cg,j,x,y):
    jx = j%(cg+1)
    jy = j//(cg+1)
    return CGbasis1d_dX(cg,jx,x)*CGbasis1d(cg,jy,y)

def CGbasisfunction_dY(cg,j,x,y):
    jx = j%(cg+1)
    jy = j//(cg+1)
    return CGbasis1d(cg,jx,x)*CGbasis1d_dX(cg,jy,y)

inversecg = [np.array([1.,1.,1.,1.]),
             np.array([9., 4.5, 9., 4.5, 2.25, 4.5, 9., 4.5, 9.0])]

