# generates c++ header files to define several static const matrices.
#
# here: projection between CG and DG vectors
#

import numpy as np


### DG
def dgdofs(d):    # Number of unknowns per element depending on gauss degree
    if d==0:
        return 1
    elif d==1:
        return 3
    elif d==2:
        return 6
    else:
        assert False,'dG3 and higher is not implemented'
    
### Inverse element mass matrix for the dg methods
inversemass = np.array([1., 12., 12., 180., 180., 144.])
        
### Gauss quadrature
gausspoints = np.array([
    [0.5,0,0,0],
    [0.5-np.sqrt(1./12.),0.5+np.sqrt(1./12.),0,0],
    [0.5-np.sqrt(3./20.),0.5,0.5+np.sqrt(3./20.),0],
    [0.5-0.5*np.sqrt(3./7.+2./7.*np.sqrt(6./5.)),
     0.5-0.5*np.sqrt(3./7.-2./7.*np.sqrt(6./5.)),
     0.5+0.5*np.sqrt(3./7.-2./7.*np.sqrt(6./5.)),
     0.5+0.5*np.sqrt(3./7.+2./7.*np.sqrt(6./5.))]])

gaussweights = np.array([
    [1.0,0,0,0],
    [0.5,0.5,0,0],
    [5./18.,8./18.,5./18.,0],
    [(18.-np.sqrt(30.0))/72.,
     (18.+np.sqrt(30.0))/72.,
     (18.+np.sqrt(30.0))/72.,
     (18.-np.sqrt(30.0))/72.]])


def sanitycheck_gauss():
    for p in range(10):    # integrate x^p for p=0,1,2,...,9
        ex = 1.0/(1.0+p)
#        print("Integrate x^{0} over [0,1]: ".format(p),end='')

        for q in range(4): # gauss 1,2,3,4
            gi = 0.0
            for k in range(q+1):
                gi = gi + gaussweights[q,k] * gausspoints[q,k]**p
            if np.fabs(gi-ex)>1.e-14:
                assert p>=2*(q+1), 'Gauss rule should be exact'

# Evaluates the dG-basis functions on [0,1]^2 in (x,y)
def DGbasisfunction(j,x,y):
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
    else:
        print("dG3 and higher not implemented (yet)")
        assert False

# Evaluates the CG(2)-basis functions on [0,1]^2 in (x,y)
def CGbasis1d(j,x):
    if j==0:
        return 2.0*(x-0.5)*(x-1.0)
    elif j==1:
        return 4.0*x*(1.0-x)
    elif j==2:
        return 2.0*x*(x-0.5)
    else:
        print("CGbasis1d only for j=0,1,2")
        assert False
# ... and its derivative
def CGbasis1d_dX(j,x):
    if j==0:
        return 4.0*x-3.0
    elif j==1:
        return 4.0-8.0*x
    elif j==2:
        return 4.0*x-1.0
    else:
        print("CGbasis1d only for j=0,1,2")
        assert False

def CGbasisfunction(j,x,y):
    jx = j%3
    jy = j//3
    return CGbasis1d(jx,x)*CGbasis1d(jy,y)
def CGbasisfunction_dX(j,x,y):
    jx = j%3
    jy = j//3
    return CGbasis1d_dX(jx,x)*CGbasis1d(jy,y)
def CGbasisfunction_dY(j,x,y):
    jx = j%3
    jy = j//3
    return CGbasis1d(jx,x)*CGbasis1d_dX(jy,y)
        
# Evaluates 1d dG-basis on the edge [0,1]
def edgebasisfunction(j,x):
    if j==0:
        return 1.
    elif j==1:
        return x-0.5
    elif j==2:
        return (x-0.5)*(x-0.5)-1.0/12.0
    else:
        print("dG3 and higher not implemented (yet)")
        assert False


# Compute the Projection matris realizing
#
# (dg, psi) = (cg, psi) in terms of dg = A * cg
#
# 'd' is the degree of the DG space
def cg2dg_matrix(d):
    # print header
    print('static const Eigen::Matrix<double, {0}, {1}, Eigen::RowMajor> CG2_to_DG{2} ='.format(dgdofs(d),9, d))
    print('\t(Eigen::Matrix<double, {0}, {1}, Eigen::RowMajor>() <<'.format(dgdofs(d),9))

    for dgi in range(dgdofs(d)):
        for cgi in range(9):
            xxx = 0
            for gx in range(3):
                for gy in range(3):
                    X = gausspoints[2][gx]
                    Y = gausspoints[2][gy]
                    xxx=xxx+gaussweights[2][gx]*gaussweights[2][gy]*CGbasisfunction(cgi,X,Y) * DGbasisfunction(dgi,X,Y)
            
            print(xxx*inversemass[dgi],end='')
            if (cgi<9-1) or dgi<dgdofs(d)-1:
                print(',',end='')
            else:
                print(').finished();')


# Compute the Projection matris realizing
#
# (dg, psi) = (d_X/Y cg, psi) in terms of dg = A_dX/Y * cg
#
# 'd' is the degree of the DG space
def cg2dg_dxy_matrix(d,dXY):
    # print header
    print('static const Eigen::Matrix<double, {0}, {1}, Eigen::RowMajor> CG2_to_DG{2}_d{3} ='.format(dgdofs(d),9, d, dXY))
    print('\t(Eigen::Matrix<double, {0}, {1}, Eigen::RowMajor>() <<'.format(dgdofs(d),9))

    for dgi in range(dgdofs(d)):
        for cgi in range(9):
            xxx = 0
            for gx in range(3):
                for gy in range(3):
                    X = gausspoints[2][gx]
                    Y = gausspoints[2][gy]
                    if (dXY=='X'):
                        xxx=xxx+gaussweights[2][gx]*gaussweights[2][gy]*CGbasisfunction_dX(cgi,X,Y) * DGbasisfunction(dgi,X,Y)
                    elif (dXY=='Y'):
                        xxx=xxx+gaussweights[2][gx]*gaussweights[2][gy]*CGbasisfunction_dY(cgi,X,Y) * DGbasisfunction(dgi,X,Y)
                    else:
                        print('Direction',dXY,'not known')
                        assert False
                        
            print(xxx*inversemass[dgi],end='')
            if (cgi<9-1) or dgi<dgdofs(d)-1:
                print(',',end='')
            else:
                print(').finished();')
                
    
#

### Main

# make sure that Guass quadrature is correct
sanitycheck_gauss()

# Some output
print('//------------------------------------')
print('#ifndef __codegeneration_cg2dg__')
print('#define __codegeneration_cg2dg__')
print('//------------------------------------')
print('\n')
print('// Automatically generated by codegeneration/project_cg_dg.py')
print('//')
print('// Generates the matrices CG2toDG[dg]')
print('// - Realizes the projection of a CG2 vector into the DG[dg] space')
print('//   dg = CG2toDG[dg] * cg')
print('')

# print out guass points and weights
print('//------------------------------ CG2toDG\n')
for dg in [1,2]:
    cg2dg_matrix(dg)
    print('')

print('\n')
print('// Generates the matrices CG2toDG[dg]_dX and CG2toDG[dg]_dY')
print('// - Realizes the projection of the derivative of a CG2 vector into the DG[dg] space')
print('//   dg = CG2toDG[dg] * cg')
print('')

# print out guass points and weights
print('//------------------------------ CG2toDG\n')
for dg in [1,2]:
    cg2dg_dxy_matrix(dg,'X')
    cg2dg_dxy_matrix(dg,'Y')
    print('')


# Some output
print('//------------------------------------')
print('#endif')
