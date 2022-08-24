# generates c++ header files to define several static const matrices
# used in the dG transport solver.
#

import numpy as np
import basisfunctions as bf
import gaussquadrature as gq

#
# evaluate basis functions in the quadrature points on  PSI
# the edges: left, right, up, down
#
# edge: left, right, up, down
# d   : number of dG - dofs (1, 3, 6, 8) for dG0/1/2 and dG2+ (nabla Q2)
# g   : number of gauss points (1,2,3)
def basisfunctions_in_gausspoints_edge(edge, d, g):
    # print header
    print('static const Eigen::Matrix<double, {0}, {1}, Eigen::RowMajor> PSI{2}{3}_{4} ='.format(g,d,d,g,edge))
    print('\t(Eigen::Matrix<double, {0}, {1}, Eigen::RowMajor>() <<'.format(g,d))
    print('\t',end=' ')
    for gp in range(g):
        for dp in range(d):
            if (edge==3):
                print(bf.inversemass[dp]*gq.gaussweights[g-1,gp] * bf.dgbasis(dp, 0.0, gq.gausspoints[g-1,gp]),end='')
            elif (edge==1):
                print(bf.inversemass[dp]*gq.gaussweights[g-1,gp] * bf.dgbasis(dp, 1.0, gq.gausspoints[g-1,gp]),end='')
            elif (edge==0):
                print(bf.inversemass[dp]*gq.gaussweights[g-1,gp] * bf.dgbasis(dp, gq.gausspoints[g-1,gp], 0.0),end='')
            elif (edge==2):
                print(bf.inversemass[dp]*gq.gaussweights[g-1,gp] * bf.dgbasis(dp, gq.gausspoints[g-1,gp], 1.0),end='')
            if (gp<g-1 or dp<d-1):
                print(', ',end='')
            else:
                print(').finished();')


# NEW! Without Mass Matrix!!!
#
# evaluate dg basis functions in the quadrature points on 
# the edges: left, right, up, down
#
# edge: left, right, up, down
# d   : number of dG - dofs (1, 3, 6, 8) for dG0/1/2 and dG2+ (nabla Q2)
# g   : number of gauss points (1,2,3)
def basisfunctions_in_gausspoints_edge_new(edge, d, g):
    # print header
    print('static const Eigen::Matrix<double, {0}, {1}, Eigen::RowMajor> PSIe_w_{2}_{3}_{4} ='.format(g,d,d,g,edge))
    print('\t(Eigen::Matrix<double, {0}, {1}, Eigen::RowMajor>() <<'.format(g,d))
    print('\t',end=' ')
    for gp in range(g):
        for dp in range(d):
            if (edge==3):
                print(gq.gaussweights[g-1,gp] * bf.dgbasis(dp, 0.0, gq.gausspoints[g-1,gp]),end='')
            elif (edge==1):
                print(gq.gaussweights[g-1,gp] * bf.dgbasis(dp, 1.0, gq.gausspoints[g-1,gp]),end='')
            elif (edge==0):
                print(gq.gaussweights[g-1,gp] * bf.dgbasis(dp, gq.gausspoints[g-1,gp], 0.0),end='')
            elif (edge==2):
                print(gq.gaussweights[g-1,gp] * bf.dgbasis(dp, gq.gausspoints[g-1,gp], 1.0),end='')
            if (gp<g-1 or dp<d-1):
                print(', ',end='')
            else:
                print(').finished();')


#
# evaluates integral with basis function in cell, scaled by inverse mass
#   M_d^(-1) ( (X), Psi_i(g) )
#      =
#   M_d^(-1) weight(g) * (X) * Psi_i(g)
#
# edge: left, right, up, down
# d   : number of dG - dofs (1, 3, 6, 8) for dG0/1/2 and dG2+ (nabla Q2)
# g   : number of gauss points (1,2,3) in each direction
def integration_basisfunctions_in_gausspoints_cell(d, g):

    
    
    # print header
    if d>1:
        print('static const Eigen::Matrix<double, {0}, {1}, Eigen::RowMajor> IBC{2}{3} ='.format(d,g*g,d,g))
        print('\t(Eigen::Matrix<double, {0}, {1}, Eigen::RowMajor>() <<'.format(d,g*g))
    else:
        print('static const Eigen::Matrix<double, {0}, {1}> IBC{2}{3} ='.format(d,g*g,d,g))
        print('\t(Eigen::Matrix<double, {0}, {1}>() <<'.format(d,g*g))
    print('\t',end=' ')
    for dp in range(d):
        for gy in range(g):
            for gx in range(g):
                print(bf.inversemass[dp]*gq.gaussweights[g-1,gx]*gq.gaussweights[g-1,gy] * bf.dgbasis(dp, gq.gausspoints[g-1,gx],gq.gausspoints[g-1,gy]),end='')
                
                if (gx<g-1 or gy<g-1 or dp<d-1):
                    print(', ',end='')
                else:
                    print(').finished();')

                
#
# evaluate basis functions in the quadrature points on
# the cell
#
# d   : number of dG - dofs (1, 3, 6, 8) for dG0/1/2 and dG2+ (nabla Q2)
# g   : Gauss points in each direction (2,3). lower left to upper right, y/x
def basisfunctions_in_gausspoints_cell(d, g):

    # print header
    print('template<> struct PSIImpl< {0}, {1} >{{'.format(d, g))
    if g>1:
        print('static inline const Eigen::Matrix<double, {0}, {1}, Eigen::RowMajor> value = (Eigen::Matrix<double, {0}, {1}, Eigen::RowMajor>() <<'.format(d,g*g))
    else:
        print('static inline const Eigen::Matrix<double, {0}, {1}> value = (Eigen::Matrix<double, {0}, {1}>() <<'.format(d,g*g))
    print('\t',end=' ')
    for dp in range(d):
        for gy in range(g):
            for gx in range(g):
                print(bf.dgbasis(dp, gq.gausspoints[g-1,gx], gq.gausspoints[g-1,gy]),end='')

                if (gx<g-1 or gy<g-1 or dp<d-1):
                    print(', ',end='')
                else:
                    print(').finished();};')



#
# evaluate basis functions in the quadrature points on
# the cell
#
# d   : number of dG - dofs (1, 3, 6, 8) for dG0/1/2 and dG2+ (nabla Q2)
# g   : Gauss points in each direction (2,3). lower left to upper right, y/x
def basisfunctions_in_gausspoints_cell_gradient(d, g):

    # print header
    print('template<> struct PSIxImpl< {0}, {1} >{{'.format(d, g))
    if g>1:
        print('static inline const Eigen::Matrix<double, {0}, {1}, Eigen::RowMajor> value = (Eigen::Matrix<double, {0}, {1}, Eigen::RowMajor>() <<'.format(d,g*g))
    else:
        print('static inline const Eigen::Matrix<double, {0}, {1}> value = (Eigen::Matrix<double, {0}, {1}>() <<'.format(d,g*g))
    print('\t',end=' ')
    for dp in range(d):
        for gy in range(g):
            for gx in range(g):
                print(bf.dx_dgbasis(dp, gq.gausspoints[g-1,gx], gq.gausspoints[g-1,gy]),end='')

                if (gx<g-1 or gy<g-1 or dp<d-1):
                    print(', ',end='')
                else:
                    print(').finished();};')

    # print header
    print('template<> struct PSIyImpl< {0}, {1} >{{'.format(d, g))
    if g>1:
        print('static inline const Eigen::Matrix<double, {0}, {1}, Eigen::RowMajor> value = (Eigen::Matrix<double, {0}, {1}, Eigen::RowMajor>() <<'.format(d,g*g))
    else:
        print('static inline const Eigen::Matrix<double, {0}, {1}> value = (Eigen::Matrix<double, {0}, {1}>() <<'.format(d,g*g))
    print('\t',end=' ')
    for dp in range(d):
        for gy in range(g):
            for gx in range(g):
                print(bf.dy_dgbasis(dp, gq.gausspoints[g-1,gx], gq.gausspoints[g-1,gy]),end='')

                if (gx<g-1 or gy<g-1 or dp<d-1):
                    print(', ',end='')
                else:
                    print(').finished();};')
                    
#
# evaluate edge basis functions in the quadrature points
# d   : number of dG - dofs (1,2,3) for dG0/1/2 and dG2+ (nabla Q2)
# g   : number of gauss points (1,2,3)
def edge_basisfunctions_in_gausspoints(d, g):

    # print header
    print('static const Eigen::Matrix<double, {0}, {1}, Eigen::RowMajor> PSIe{2}{3} ='.format(d,g,d,g))
    print('\t(Eigen::Matrix<double, {0}, {1}, Eigen::RowMajor>() <<'.format(d,g))
    print('\t',end=' ')
    for dp in range(d):
        for gp in range(g):
            print(bf.dgbasis_edge(dp, gq.gausspoints[g-1,gp]),end='')
            if (dp<d-1 or gp<g-1):
                print(', ',end='')
            else:
                print(').finished();')





### Main

# make sure that Guass quadrature is correct
gq.sanitycheck_gauss()

# Some output
print('#ifndef __BASISFUNCTIONSGUASSPOINTS_HPP')
print('#define __BASISFUNCTIONSGUASSPOINTS_HPP')
print('\n')
print('// Automatically generated by codegeneration/basisfunctions_gausspoints.py')
print('//')
print('// Generates the vectors gauss_points[gq] and gauss_weights[gq]')
print('// - stores the points and weights of the gq-point Guass rule')
print('// - the integration is scaled to [0,1]')
print('')
print('// Generates the matrices PSI[dg][gq]_[e]')
print('// - dg is the degree of the dG space')
print('// - gq is the number of Gauss quadrature points')
print('// - e is the edge number (0-lower, 1-right, 2-up, 3-left)')
print('//')
print('// - Each matrix PSI_e[i,j] stores the value of the j-th basis function')
print('//   in the i-th Gauss quadrature point alongt the edge e, weighted')
print('//   with corresponding Gauss weight and the inverse of the mass matrix')
print('//   bf.inversemass(j) * phi_j( gauss_point(i) ) * gauss_weight(i)')
print('')
print('// Generates the matrices PSIe[dg][gq]')
print('// - stores the value of the basis functions on the edge in ')
print('//   the Guass points along the edge')
print('// - dg is the dG degree and gq the number of Gauss points')
print('// - PSIe[i,j] is simply phi_j( gauss_point(i) )')
print('')
print('// Generates the matrices PSI[dg][gq]')
print('// - stores the value of the basis functions on the cell in ')
print('//   the Guass points')
print('// - dg is the dG degree and gq the number of Gauss points in each direction')
print('// - PSI[i,j] is phi_j( gauss_point(i) )')
print('')


# print out guass points and weights
print('\n\n//------------------------------ Gauss Quadrature\n')
for gp in [1,2,3]:
    print('constexpr double gauss_points{0}[{0}] = {{'.format(gp),end='')
    for q in range(gp):
        print(gq.gausspoints[gp-1,q],end='')
        if q<gp-1:
            print(',',end='')
    print('};')
    print('constexpr double gauss_weights{0}[{0}] = {{'.format(gp),end='')
    for q in range(gp):
        print(gq.gaussweights[gp-1,q],end='')
        if q<gp-1:
            print(',',end='')
    print('};')
    
print('\n\n//------------------------------ Gauss Quadrature\n')
for g in [1,2,3]:
    print('static const Eigen::Matrix<double, 1, {0}, Eigen::RowMajor> GAUSSWEIGHTS_{1} ='.format(g*g,g))
    print('\t(Eigen::Matrix<double, 1, {0}, Eigen::RowMajor>() <<'.format(g*g))
    print('\t',end=' ')
    for gy in range(g):
        for gx in range(g):
            print(gq.gaussweights[g-1][gx]*gq.gaussweights[g-1][gy],end='')
            
            if (gx<g-1 or gy<g-1):
                print(', ',end='')
            else:
                print(').finished();')


print('\n\n//------------------------------ Gauss Points\n')
for g in [1,2,3]:
    if g>1:
        print('static const Eigen::Matrix<double, 2, {0}, Eigen::RowMajor> GAUSSPOINTS_{1} ='.format(g*g,g))
        print('\t(Eigen::Matrix<double, 2, {0}, Eigen::RowMajor>() <<'.format(g*g))
    else:
        print('static const Eigen::Matrix<double, 2, {0}> GAUSS_{1} ='.format(g*g,g))
        print('\t(Eigen::Matrix<double, 2, {0}>() <<'.format(g*g))
    print('\t',end=' ')
    for gy in range(g):
        for gx in range(g):
            print(gq.gausspoints[g-1][gx],end=',')
    for gy in range(g):
        for gx in range(g):
            print(gq.gausspoints[g-1][gy],end='')
            
            if (gx<g-1 or gy<g-1):
                print(', ',end='')
            else:
                print(').finished();')
    


print('\n\n//------------------------------ Basis Functions in Gauss Points (edge)\n')
# generate arrays that evaluate basis functions in the gauss points on the edges
for dg in [1,3,6,8]:
    for e in [0,1,2,3]:
        if dg==1:
            basisfunctions_in_gausspoints_edge(e, dg,1)
        elif dg==3:
            basisfunctions_in_gausspoints_edge(e, dg,2)
        else:
            basisfunctions_in_gausspoints_edge(e, dg,3)
        print('')

        if dg==1:
            basisfunctions_in_gausspoints_edge_new(e, dg,1)
        elif dg==3:
            basisfunctions_in_gausspoints_edge_new(e, dg,2)
        else:
            basisfunctions_in_gausspoints_edge_new(e, dg,3)
        print('')
        

print('\n\n//------------------------------ Basis Functions in Gauss Points (cell)\n')

print('template<int DG, int GP>')
print('struct PSIImpl;')
for dg in [1, 3,6,8]:
    basisfunctions_in_gausspoints_cell(dg,1)
    basisfunctions_in_gausspoints_cell(dg,2)
    basisfunctions_in_gausspoints_cell(dg,3)
print('template<int DG, int GP>')
print('const Eigen::Matrix<double, DG, GP*GP, Eigen::RowMajor> PSI = PSIImpl<DG, GP>::value;')
print('')


    
print('template<int DG, int GP>')
print('struct PSIxImpl;')
print('template<int DG, int GP>')
print('struct PSIyImpl;')    
for dg in [1, 3,6,8]:
    basisfunctions_in_gausspoints_cell_gradient(dg,1)
    basisfunctions_in_gausspoints_cell_gradient(dg,2)
    basisfunctions_in_gausspoints_cell_gradient(dg,3)
print('template<int DG, int GP>')
print('const Eigen::Matrix<double, DG, GP*GP, Eigen::RowMajor> PSIx = PSIxImpl<DG, GP>::value;')
print('template<int DG, int GP>')
print('const Eigen::Matrix<double, DG, GP*GP, Eigen::RowMajor> PSIy = PSIyImpl<DG, GP>::value;')
print('')

    

for dg in [1,3,6,8]:
    integration_basisfunctions_in_gausspoints_cell(dg,2)
    integration_basisfunctions_in_gausspoints_cell(dg,3)

    
print('\n\n//------------------------------ Edge Basis Functions in Gauss Points\n')
# generate arrays that evaluate basis functions in the gauss points on the edges
for dg in [1,2,3]:
    edge_basisfunctions_in_gausspoints(dg,dg)
    print('')

edge_basisfunctions_in_gausspoints(2,3)
print('') 

# Some output
print('#endif /* __BASISFUNCTIONSGUASSPOINTS_HPP */')
