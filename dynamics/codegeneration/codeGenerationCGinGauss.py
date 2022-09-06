# generates c++ header files to define several static const matrices
# used in the dG transport solver.
#

import numpy as np
import basisfunctions as bf
import gaussquadrature as gq
                




# evaluate cg basis functions in the quadrature points on
# the cell
#
# cg  : cg degree (1 or 2)
# gp  : Gauss points in each direction (2,3). lower left to upper right, y/x
#
# PHI
def cgbasisfunctions_in_gausspoints(cg, gp):

    # print header
    print('template<> struct PHIImpl< {0}, {1} >{{'.format(cg, gp))
    if gp>1:
        print('static inline const Eigen::Matrix<double, {0}, {1}, Eigen::RowMajor> value = (Eigen::Matrix<double, {0}, {1}, Eigen::RowMajor>() <<'.format(bf.cgdofs(cg),gp*gp))
    else:
        print('static inline const Eigen::Matrix<double, {0}, {1}> value = (Eigen::Matrix<double, {0}, {1}>() <<'.format(bf.cgdofs(cg),gp*gp))
    print('\t',end=' ')

    for dp in range(bf.cgdofs(cg)):
        for gy in range(gp):
            for gx in range(gp):
                print(bf.CGbasisfunction(cg, dp, gq.gausspoints[gp-1,gx], gq.gausspoints[gp-1,gy]),end='')

                if (gx<gp-1 or gy<gp-1 or dp<bf.cgdofs(cg)-1):
                    print(', ',end='')
                else:
                    print(').finished();};')

                    

def cg_dx_basisfunctions_in_gausspoints(cg, gp):
    # print header
    print('template<> struct PHIxImpl< {0}, {1} >{{'.format(cg, gp))
    if gp>1:
        print('static inline const Eigen::Matrix<double, {0}, {1}, Eigen::RowMajor> value = (Eigen::Matrix<double, {0}, {1}, Eigen::RowMajor>() <<'.format(bf.cgdofs(cg),gp*gp))
    else:
        print('static inline const Eigen::Matrix<double, {0}, {1}> value = (Eigen::Matrix<double, {0}, {1}>() <<'.format(bf.cgdofs(cg),gp*gp))
    print('\t',end=' ')
    for dp in range(bf.cgdofs(cg)):
        for gy in range(gp):
            for gx in range(gp):
                print(bf.CGbasisfunction_dX(cg, dp, gq.gausspoints[gp-1,gx], gq.gausspoints[gp-1,gy]),end='')

                if (gx<gp-1 or gy<gp-1 or dp<bf.cgdofs(cg)-1):
                    print(', ',end='')
                else:
                    print(').finished();};')

def cg_dy_basisfunctions_in_gausspoints(cg, gp):
    # print header
    print('template<> struct PHIyImpl< {0}, {1} >{{'.format(cg, gp))
    if gp>1:
        print('static inline const Eigen::Matrix<double, {0}, {1}, Eigen::RowMajor> value = (Eigen::Matrix<double, {0}, {1}, Eigen::RowMajor>() <<'.format(bf.cgdofs(cg),gp*gp))
    else:
        print('static inline const Eigen::Matrix<double, {0}, {1}> value = (Eigen::Matrix<double, {0}, {1}>() <<'.format(bf.cgdofs(cg),gp*gp))
    print('\t',end=' ')

    for dp in range(bf.cgdofs(cg)):
        for gy in range(gp):
            for gx in range(gp):
                print(bf.CGbasisfunction_dY(cg, dp, gq.gausspoints[gp-1,gx], gq.gausspoints[gp-1,gy]),end='')

                if (gx<gp-1 or gy<gp-1 or dp<bf.cgdofs(cg)-1):
                    print(', ',end='')
                else:
                    print(').finished();};')


# evaluate a cg function in the quadrature points
# cg  : cg degree (1 or 2)
# gp  : Gauss points in each direction (2,3). lower left to upper right, y/x
# matrix of size cgdofs x gausspoints (2d)
def cgfunction_in_gausspoints(cg, gp):

    # print header
    if gp>1:
        print('static const Eigen::Matrix<double, {0}, {1}, Eigen::RowMajor> CG_CG{2}FUNC_in_GAUSS{3} ='.format(gp*gp,bf.cgdofs(cg),cg,gp))
        print('\t(Eigen::Matrix<double, {0}, {1}, Eigen::RowMajor>() <<'.format(gp*gp,bf.cgdofs(cg)))
    else:
        print('static const Eigen::Matrix<double, {0}, {1}> CG_CG{2}FUNC_in_GAUSS{3} ='.format(gp*gp,bf.cgdofs(cg), cg, gp))
        print('\t(Eigen::Matrix<double, {0}, {1}>() <<'.format(gp*gp,bf.cgdofs(cg)))
    print('\t',end=' ')
    for gy in range(gp):
        for gx in range(gp):
            for dp in range(bf.cgdofs(cg)):

                print(bf.CGbasisfunction(cg, dp, gq.gausspoints[gp-1,gx], gq.gausspoints[gp-1,gy]),end='')

                if (gx<gp-1 or gy<gp-1 or dp<bf.cgdofs(cg)-1):
                    print(', ',end='')
                else:
                    print(').finished();')


# evaluate the derivative of a cg function in the quadrature points
# cg  : cg degree (1 or 2)
# gp  : Gauss points in each direction (2,3). lower left to upper right, y/x
# matrix of size cgdofs x gausspoints (2d)
def cgfunction_dx_in_gausspoints(cg, gp):

    # print header
    if gp>1:
        print('static const Eigen::Matrix<double, {0}, {1}, Eigen::RowMajor> CG_CG{2}FUNC_DX_in_GAUSS{3} ='.format(gp*gp,bf.cgdofs(cg),cg,gp))
        print('\t(Eigen::Matrix<double, {0}, {1}, Eigen::RowMajor>() <<'.format(gp*gp,bf.cgdofs(cg)))
    else:
        print('static const Eigen::Matrix<double, {0}, {1}> CG_CG{2}FUNC_DX_in_GAUSS{3} ='.format(gp*gp,bf.cgdofs(cg), cg, gp))
        print('\t(Eigen::Matrix<double, {0}, {1}>() <<'.format(gp*gp,bf.cgdofs(cg)))
    print('\t',end=' ')
    for gy in range(gp):
        for gx in range(gp):
            for dp in range(bf.cgdofs(cg)):

                print(bf.CGbasisfunction_dX(cg, dp, gq.gausspoints[gp-1,gx], gq.gausspoints[gp-1,gy]),end='')

                if (gx<gp-1 or gy<gp-1 or dp<bf.cgdofs(cg)-1):
                    print(', ',end='')
                else:
                    print(').finished();')

# evaluate the derivative of a cg function in the quadrature points
# cg  : cg degree (1 or 2)
# gp  : Gauss points in each direction (2,3). lower left to upper right, y/x
# matrix of size cgdofs x gausspoints (2d)
def cgfunction_dy_in_gausspoints(cg, gp):

    # print header
    if gp>1:
        print('static const Eigen::Matrix<double, {0}, {1}, Eigen::RowMajor> CG_CG{2}FUNC_DY_in_GAUSS{3} ='.format(gp*gp,bf.cgdofs(cg),cg,gp))
        print('\t(Eigen::Matrix<double, {0}, {1}, Eigen::RowMajor>() <<'.format(gp*gp,bf.cgdofs(cg)))
    else:
        print('static const Eigen::Matrix<double, {0}, {1}> CG_CG{2}FUNC_DY_in_GAUSS{3} ='.format(gp*gp,bf.cgdofs(cg), cg, gp))
        print('\t(Eigen::Matrix<double, {0}, {1}>() <<'.format(gp*gp,bf.cgdofs(cg)))
    print('\t',end=' ')
    for gy in range(gp):
        for gx in range(gp):
            for dp in range(bf.cgdofs(cg)):

                print(bf.CGbasisfunction_dY(cg, dp, gq.gausspoints[gp-1,gx], gq.gausspoints[gp-1,gy]),end='')

                if (gx<gp-1 or gy<gp-1 or dp<bf.cgdofs(cg)-1):
                    print(', ',end='')
                else:
                    print(').finished();')



### Main

# make sure that Guass quadrature is correct
gq.sanitycheck_gauss()

# Some output
print('#ifndef __CODEGENERATION_CG_IN_GAUSS_HPP')
print('#define __CODEGENERATION_CG_IN_GAUSS_HPP')
print('\n')
print('// Automatically generated by codegeneration/codeGenerationCGinGauss.py')
print('//')
print('// Generates the matrices CG_CG[1/2]_in_GAUSS[1/2/3] ')
print('// - matrices of size ndofs x ngausspoints')
print('// - X_{iq} = phi_i(xq) ')
print('')
print('')

print('// size of local cgbasis')
print('#define CGDOFS(CG) ( (CG==1)?4:9 )')

print('\n\n//------------------------------ CG_CG[1/2]_in_GAUSS[1/2/3]\n')
# generate arrays that evaluate basis functions in the gauss points


print('template<int CG, int GP>')
print('struct PHIImpl;')
print('template<int CG, int GP>')
print('struct PHIxImpl;')
print('template<int CG, int GP>')
print('struct PHIyImpl;')
for cg in [1,2]:
    for gp in [1,2,3]:
        cgbasisfunctions_in_gausspoints(cg,gp)

        cg_dx_basisfunctions_in_gausspoints(cg,gp)

        cg_dy_basisfunctions_in_gausspoints(cg,gp)
print('template<int CG, int GP>')
print('const Eigen::Matrix<double, CGDOFS(CG), GP*GP, Eigen::RowMajor> PHI = PHIImpl<CG, GP>::value;')
print('template<int CG, int GP>')
print('const Eigen::Matrix<double, CGDOFS(CG), GP*GP, Eigen::RowMajor> PHIx = PHIxImpl<CG, GP>::value;')
print('template<int CG, int GP>')
print('const Eigen::Matrix<double, CGDOFS(CG), GP*GP, Eigen::RowMajor> PHIy = PHIyImpl<CG, GP>::value;')
print('')



print('')
print('\n\n//------------------------------ CG_CG[1/2]FUNC_in_GAUSS[1/2/3]\n')
# generate arrays that evaluates a cg function in the gauss points
for cg in [1,2]:
    for gp in [1,2,3]:
        cgfunction_in_gausspoints(cg,gp)


print('')
print('\n\n//------------------------------ CG_CG[1/2]FUNC_in_GAUSS[1/2/3]\n')
# generate arrays that evaluates the derivatives of a cg functions in the gauss points
for cg in [1,2]:
    for gp in [1,2,3]:
        cgfunction_dx_in_gausspoints(cg,gp)
        cgfunction_dy_in_gausspoints(cg,gp)

# Some output
print('#endif /* __BASISFUNCTIONSGUASSPOINTS_HPP */')
