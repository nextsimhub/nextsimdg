# defines the dg and cg basis functions on the unit square [0,1]^2
# used in the different code generation functions
#

import numpy as np


# returns the number of local unknowns in the dg spaces
def dgdofs(d):    # Number of unknowns per element depending on gauss degree
    if d==0:
        return 1
    elif d==1:
        return 3
    elif d==2:
        return 6
    else:
        assert False,'dG3 and higher is not implemented'

# returns the number of local unknowns in the cg spaces
def cgdofs(d):    # Number of unknowns per element depending on gauss degree
    if d==1:
        return 4
    elif d==2:
        return 9
    else:
        assert False,'only cg1 and cg2 are supported'

# Evaluates the dG-basis functions on [0,1]^2 in (x,y)
# dG 0:  1
# dG 1:  x, y
# dG 2:  x^2, y^2, xy
# dG X:  x^2y, xy^2
def dgbasis(j,x,y):
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
        assert False

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
        assert False

def dy_dgbasis(j,x,y):
    if j==0:
        return 0.
    elif j==1:
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
        assert False


        
# Evaluates 1d dG-basis on the edge [0,1]
def dgbasis_edge(j,x):
    if j==0:
        return 1.
    elif j==1:
        return x-0.5
    elif j==2:
        return (x-0.5)*(x-0.5)-1.0/12.0
    else:
        print("dG3 and higher not implemented (yet)")
        assert False

### Inverse element mass matrix for the dg methods
inversemass = np.array([1., 12., 12., 180., 180., 144., 2160., 2160.])
        

# Evaluates the CG(2)-basis functions on [0,1]^2 in (x,y)
def CGbasis1d(cg,j,x):
    if cg==1:
        if j==0:
            return 1.0-x
        elif j==1:
            return x
        else:
            print("CG1basis1d only for j=0,1")
            assert False

    elif cg==2:
        if j==0:
            return 2.0*(x-0.5)*(x-1.0)
        elif j==1:
            return 4.0*x*(1.0-x)
        elif j==2:
            return 2.0*x*(x-0.5)
        else:
            print("CGbasis1d only for j=0,1,2")
            assert False
    else:
        print("only CG1 CG2")
        assert False

# ... and its derivative
def CGbasis1d_dX(cg,j,x):
    if cg==1:
        if j==0:
            return -1
        else:
            return 1
    else:
        if j==0:
            return 4.0*x-3.0
        elif j==1:
            return 4.0-8.0*x
        elif j==2:
            return 4.0*x-1.0
        else:
            print("CGbasis1d only for j=0,1,2")
            assert False

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

