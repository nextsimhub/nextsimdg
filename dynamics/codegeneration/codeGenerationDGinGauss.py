# generates c++ header files to define several static const matrices
# used in the dG transport solver.
#

import basisfunctions as bf
import gaussquadrature as gq
import numpy as np


#
# evaluate basis functions in the quadrature points on  PSI
# the edges: left, right, up, down
#
# edge: left, right, up, down
# d   : number of dG - dofs (1, 3, 6, 8) for dG0/1/2 and dG2+ (nabla Q2)
# g   : number of gauss points (1,2,3)
# dtype: C++ type name ("float" or "double")
def basisfunctions_in_gausspoints_edge(edge, d, g, dtype):
    print(
        f"static const Eigen::Matrix<{dtype}, {g}, {d}, Eigen::RowMajor> PSI{d}{g}{dtype}_{edge} ="
    )
    if g > 1 and d > 1:
        print(f"\t(Eigen::Matrix<{dtype}, {g}, {d}, Eigen::RowMajor>() <<")
    else:
        print(f"\t(Eigen::Matrix<{dtype}, {g}, {d}>() <<")

    print("\t", end=" ")
    for gp in range(g):
        for dp in range(d):
            if edge == 3:
                print(
                    bf.inversemass[dp]
                    * gq.gaussweights[g - 1, gp]
                    * bf.dgbasis(dp, 0.0, gq.gausspoints[g - 1, gp]),
                    end="",
                )
            elif edge == 1:
                print(
                    bf.inversemass[dp]
                    * gq.gaussweights[g - 1, gp]
                    * bf.dgbasis(dp, 1.0, gq.gausspoints[g - 1, gp]),
                    end="",
                )
            elif edge == 0:
                print(
                    bf.inversemass[dp]
                    * gq.gaussweights[g - 1, gp]
                    * bf.dgbasis(dp, gq.gausspoints[g - 1, gp], 0.0),
                    end="",
                )
            elif edge == 2:
                print(
                    bf.inversemass[dp]
                    * gq.gaussweights[g - 1, gp]
                    * bf.dgbasis(dp, gq.gausspoints[g - 1, gp], 1.0),
                    end="",
                )
            if gp < g - 1 or dp < d - 1:
                print(", ", end="")
            else:
                print(").finished();")


# NEW! Without Mass Matrix!!!
#
# evaluate dg basis functions in the quadrature points on
# the edges: left, right, up, down
#
# edge: left, right, up, down
# d   : number of dG - dofs (1, 3, 6, 8) for dG0/1/2 and dG2+ (nabla Q2)
# g   : number of gauss points (1,2,3)
# dtype: C++ type name ("float" or "double")
def basisfunctions_in_gausspoints_edge_new(edge, d, g, dtype):
    if g > 1 and d > 1:
        print(f"template<> struct PSIe_wImpl< {d}, {g}, {edge}, {dtype} >{{")
        print(
            f"static inline const Eigen::Matrix<{dtype}, {g}, {d}, Eigen::RowMajor> value = (Eigen::Matrix<{dtype}, {g}, {d}, Eigen::RowMajor>() <<"
        )
    else:
        print(f"template<> struct PSIe_wImpl< {d}, {g}, {edge}, {dtype} >{{")
        print(
            f"static inline const Eigen::Matrix<{dtype}, {g}, {d}> value = (Eigen::Matrix<{dtype}, {g}, {d}>() <<"
        )
    print("\t", end=" ")

    for gp in range(g):
        for dp in range(d):
            if edge == 3:
                print(
                    gq.gaussweights[g - 1, gp]
                    * bf.dgbasis(dp, 0.0, gq.gausspoints[g - 1, gp]),
                    end="",
                )
            elif edge == 1:
                print(
                    gq.gaussweights[g - 1, gp]
                    * bf.dgbasis(dp, 1.0, gq.gausspoints[g - 1, gp]),
                    end="",
                )
            elif edge == 0:
                print(
                    gq.gaussweights[g - 1, gp]
                    * bf.dgbasis(dp, gq.gausspoints[g - 1, gp], 0.0),
                    end="",
                )
            elif edge == 2:
                print(
                    gq.gaussweights[g - 1, gp]
                    * bf.dgbasis(dp, gq.gausspoints[g - 1, gp], 1.0),
                    end="",
                )
            if gp < g - 1 or dp < d - 1:
                print(", ", end="")
            else:
                print(").finished();};")


#
# evaluates integral with basis function in cell, scaled by inverse mass
#   M_d^(-1) ( (X), Psi_i(g) )
#      =
#   M_d^(-1) weight(g) * (X) * Psi_i(g)
#
# edge: left, right, up, down
# d   : number of dG - dofs (1, 3, 6, 8) for dG0/1/2 and dG2+ (nabla Q2)
# g   : number of gauss points (1,2,3) in each direction
# dtype: C++ type name ("float" or "double")
def integration_basisfunctions_in_gausspoints_cell(d, g, dtype):
    if d > 1:
        print(
            f"static const Eigen::Matrix<{dtype}, {d}, {g * g}, Eigen::RowMajor> IBC{d}{g}{dtype} ="
        )
        print(f"\t(Eigen::Matrix<{dtype}, {d}, {g * g}, Eigen::RowMajor>() <<")
    else:
        print(f"static const Eigen::Matrix<{dtype}, {d}, {g * g}> IBC{d}{g}{dtype} =")
        print(f"\t(Eigen::Matrix<{dtype}, {d}, {g * g}>() <<")
    print("\t", end=" ")
    for dp in range(d):
        for gy in range(g):
            for gx in range(g):
                print(
                    bf.inversemass[dp]
                    * gq.gaussweights[g - 1, gx]
                    * gq.gaussweights[g - 1, gy]
                    * bf.dgbasis(
                        dp, gq.gausspoints[g - 1, gx], gq.gausspoints[g - 1, gy]
                    ),
                    end="",
                )

                if gx < g - 1 or gy < g - 1 or dp < d - 1:
                    print(", ", end="")
                else:
                    print(").finished();")


#
# evaluate basis functions in the trapecoidal / simpson points on
# the cell
#
# d   : number of dG - dofs (1, 3, 6, 8) for dG0/1/2 and dG2+ (nabla Q2)
# g   : 2 Trapecoidal, 3 Simpson
# dtype: C++ type name ("float" or "double")
def basisfunctions_in_lagrangepoints_cell(d, g, dtype):
    print(f"template<> struct PSILagrangeImpl< {d}, {g}, {dtype} >{{")
    print(
        f"static inline const Eigen::Matrix<{dtype}, {d}, {g * g}, Eigen::RowMajor> value = (Eigen::Matrix<{dtype}, {d}, {g * g}, Eigen::RowMajor>() <<"
    )
    print("\t", end=" ")
    for dp in range(d):
        for gy in range(g):
            for gx in range(g):
                print(
                    bf.dgbasis(
                        dp, gq.lagrangepoints[g - 2, gx], gq.lagrangepoints[g - 2, gy]
                    ),
                    end="",
                )

                if gx < g - 1 or gy < g - 1 or dp < d - 1:
                    print(", ", end="")
                else:
                    print(").finished();};")


#
# evaluate basis functions in the quadrature points on
# the cell
#
# d   : number of dG - dofs (1, 3, 6, 8) for dG0/1/2 and dG2+ (nabla Q2)
# g   : Gauss points in each direction (2,3). lower left to upper right, y/x
# dtype: C++ type name ("float" or "double")
def basisfunctions_in_gausspoints_cell(d, g, dtype):
    print(f"template<> struct PSIImpl< {d}, {g}, {dtype} >{{")
    if g > 1:
        print(
            f"static inline const Eigen::Matrix<{dtype}, {d}, {g * g}, Eigen::RowMajor> value = (Eigen::Matrix<{dtype}, {d}, {g * g}, Eigen::RowMajor>() <<"
        )
    else:
        print(
            f"static inline const Eigen::Matrix<{dtype}, {d}, {g * g}> value = (Eigen::Matrix<{dtype}, {d}, {g * g}>() <<"
        )
    print("\t", end=" ")
    for dp in range(d):
        for gy in range(g):
            for gx in range(g):
                print(
                    bf.dgbasis(
                        dp, gq.gausspoints[g - 1, gx], gq.gausspoints[g - 1, gy]
                    ),
                    end="",
                )

                if gx < g - 1 or gy < g - 1 or dp < d - 1:
                    print(", ", end="")
                else:
                    print(").finished();};")


#
# evaluate basis functions in the quadrature points on
# the cell
#
# d   : number of dG - dofs (1, 3, 6, 8) for dG0/1/2 and dG2+ (nabla Q2)
# g   : Gauss points in each direction (2,3). lower left to upper right, y/x
# dtype: C++ type name ("float" or "double")
def basisfunctions_in_gausspoints_cell_gradient(d, g, dtype):
    print(f"template<> struct PSIxImpl< {d}, {g}, {dtype} >{{")
    if g > 1:
        print(
            f"static inline const Eigen::Matrix<{dtype}, {d}, {g * g}, Eigen::RowMajor> value = (Eigen::Matrix<{dtype}, {d}, {g * g}, Eigen::RowMajor>() <<"
        )
    else:
        print(
            f"static inline const Eigen::Matrix<{dtype}, {d}, {g * g}> value = (Eigen::Matrix<{dtype}, {d}, {g * g}>() <<"
        )
    print("\t", end=" ")
    for dp in range(d):
        for gy in range(g):
            for gx in range(g):
                print(
                    bf.dx_dgbasis(
                        dp, gq.gausspoints[g - 1, gx], gq.gausspoints[g - 1, gy]
                    ),
                    end="",
                )

                if gx < g - 1 or gy < g - 1 or dp < d - 1:
                    print(", ", end="")
                else:
                    print(").finished();};")

    print(f"template<> struct PSIyImpl< {d}, {g}, {dtype} >{{")
    if g > 1:
        print(
            f"static inline const Eigen::Matrix<{dtype}, {d}, {g * g}, Eigen::RowMajor> value = (Eigen::Matrix<{dtype}, {d}, {g * g}, Eigen::RowMajor>() <<"
        )
    else:
        print(
            f"static inline const Eigen::Matrix<{dtype}, {d}, {g * g}> value = (Eigen::Matrix<{dtype}, {d}, {g * g}>() <<"
        )
    print("\t", end=" ")
    for dp in range(d):
        for gy in range(g):
            for gx in range(g):
                print(
                    bf.dy_dgbasis(
                        dp, gq.gausspoints[g - 1, gx], gq.gausspoints[g - 1, gy]
                    ),
                    end="",
                )

                if gx < g - 1 or gy < g - 1 or dp < d - 1:
                    print(", ", end="")
                else:
                    print(").finished();};")


#
# evaluate edge basis functions in the quadrature points
# d   : number of dG - dofs (1,2,3) for dG0/1/2 and dG2+ (nabla Q2)
# g   : number of gauss points (1,2,3)
# dtype: C++ type name ("float" or "double")
def edge_basisfunctions_in_gausspoints(d, g, dtype):
    print(f"template<> struct PSIeImpl< {d}, {g}, {dtype} >{{")
    if g > 1 and d > 1:
        print(
            f"static inline const Eigen::Matrix<{dtype}, {d}, {g}, Eigen::RowMajor> value = (Eigen::Matrix<{dtype}, {d}, {g}, Eigen::RowMajor>() <<"
        )
    else:
        print(
            f"static inline const Eigen::Matrix<{dtype}, {d}, {g}> value = (Eigen::Matrix<{dtype}, {d}, {g}>() <<"
        )
    print("\t", end=" ")
    for dp in range(d):
        for gp in range(g):
            print(bf.dgbasis_edge(dp, gq.gausspoints[g - 1, gp]), end="")
            if dp < d - 1 or gp < g - 1:
                print(", ", end="")
            else:
                print(").finished();};")


### Main

# make sure that Guass quadrature is correct
gq.sanitycheck_gauss()

include_guard = "__CODEGENERATIONDGINGAUSS_HPP"

# Some output
print(f"#ifndef {include_guard}")
print(f"#define {include_guard}")
print("\n")
print("// Automatically generated by codegeneration/basisfunctions_gausspoints.py")
print("//")
print("// Generates the vectors gauss_points[gq] and gauss_weights[gq]")
print("// - stores the points and weights of the gq-point Guass rule")
print("// - the integration is scaled to [0,1]")
print("")
print("// Generates the matrices PSI[dg][gq]_[e]")
print("// - dg is the degree of the dG space")
print("// - gq is the number of Gauss quadrature points")
print("// - e is the edge number (0-lower, 1-right, 2-up, 3-left)")
print("//")
print("// - Each matrix PSI_e[i,j] stores the value of the j-th basis function")
print("//   in the i-th Gauss quadrature point alongt the edge e, weighted")
print("//   with corresponding Gauss weight and the inverse of the mass matrix")
print("//   bf.inversemass(j) * phi_j( gauss_point(i) ) * gauss_weight(i)")
print("")
print("// Generates the matrices PSIe[dg][gq]")
print("// - stores the value of the basis functions on the edge in ")
print("//   the Guass points along the edge")
print("// - dg is the dG degree and gq the number of Gauss points")
print("// - PSIe[i,j] is simply phi_j( gauss_point(i) )")
print("")
print("// Generates the matrices PSI[dg][gq]")
print("// - stores the value of the basis functions on the cell in ")
print("//   the Guass points")
print("// - dg is the dG degree and gq the number of Gauss points in each direction")
print("// - PSI[i,j] is phi_j( gauss_point(i) )")
print("")

print('#include "NextsimDynamics.hpp"')
print("\nnamespace Nextsim {")

# print out guass points and weights
print("\n\n//------------------------------ Gauss Quadrature\n")
for gp in [1, 2, 3, 4]:
    print("constexpr FloatType gauss_points{0}[{0}] = {{".format(gp), end="")
    for q in range(gp):
        print(gq.gausspoints[gp - 1, q], end="")
        if q < gp - 1:
            print(",", end="")
    print("};")
    print("constexpr FloatType gauss_weights{0}[{0}] = {{".format(gp), end="")
    for q in range(gp):
        print(gq.gaussweights[gp - 1, q], end="")
        if q < gp - 1:
            print(",", end="")
    print("};")

print("\n\n//------------------------------ Gauss Quadrature\n")

print("template<int NGP, typename T = FloatType>")
print("struct GAUSSWEIGHTSImpl;")
for g in [1, 2, 3, 4]:
    for dtype in ["float", "double"]:
        print("template<> struct GAUSSWEIGHTSImpl< {0}, {1} >{{".format(g, dtype))
        print(
            f"static inline const Eigen::Matrix<{dtype}, 1, {g * g}, Eigen::RowMajor> value = (Eigen::Matrix<{dtype}, 1, {g * g}, Eigen::RowMajor>() <<"
        )
        print("\t", end=" ")

        for gy in range(g):
            for gx in range(g):
                print(gq.gaussweights[g - 1][gx] * gq.gaussweights[g - 1][gy], end="")

                if gx < g - 1 or gy < g - 1:
                    print(", ", end="")
                else:
                    print(").finished();};")
print("template<int NGP, typename T = FloatType>")
print(
    "const Eigen::Matrix<T, 1, NGP*NGP, Eigen::RowMajor> GAUSSWEIGHTS = GAUSSWEIGHTSImpl<NGP, T>::value;"
)
print("")


print("\n\n//------------------------------ Gauss Points\n")
for g in [1, 2, 3, 4]:
    if g > 1:
        print(
            f"static const Eigen::Matrix<FloatType, 2, {g * g}, Eigen::RowMajor> GAUSSPOINTS_{g} ="
        )
        print(f"\t(Eigen::Matrix<FloatType, 2, {g * g}, Eigen::RowMajor>() <<")
    else:
        print(f"static const Eigen::Matrix<FloatType, 2, {g * g}> GAUSS_{g} =")
        print(f"\t(Eigen::Matrix<FloatType, 2, {g * g}>() <<")
    print("\t", end=" ")
    for _ in range(g):
        for gx in range(g):
            print(gq.gausspoints[g - 1][gx], end=",")
    for gy in range(g):
        for gx in range(g):
            print(gq.gausspoints[g - 1][gy], end="")

            if gx < g - 1 or gy < g - 1:
                print(", ", end="")
            else:
                print(").finished();")


print("\n\n//------------------------------ Basis Functions in Gauss Points (edge)\n")
print("template<int DG, int GP, int E, typename T = FloatType>")
print("struct PSIe_wImpl;")
for dg in [1, 3, 6, 8]:
    for e in [0, 1, 2, 3]:
        if dg == 1:
            for dtype in ["float", "double"]:
                basisfunctions_in_gausspoints_edge(e, dg, 1, dtype)
        elif dg == 3:
            for dtype in ["float", "double"]:
                basisfunctions_in_gausspoints_edge(e, dg, 2, dtype)
        else:
            for dtype in ["float", "double"]:
                basisfunctions_in_gausspoints_edge(e, dg, 3, dtype)
        print("")

        if dg == 1:
            for dtype in ["float", "double"]:
                basisfunctions_in_gausspoints_edge_new(e, dg, 1, dtype)
        elif dg == 3:
            for dtype in ["float", "double"]:
                basisfunctions_in_gausspoints_edge_new(e, dg, 2, dtype)
        else:
            for dtype in ["float", "double"]:
                basisfunctions_in_gausspoints_edge_new(e, dg, 3, dtype)
        print("")
print("template<int DG, int GP, int E, typename T = FloatType>")
print(
    "const Eigen::Matrix<T, GP, DG, Eigen::RowMajor> PSIe_w = PSIe_wImpl<DG, GP, E, T>::value;"
)
print("")


print(
    "\n\n//------------------------------ Basis Functions in Lagrange Points (cell)\n"
)

print("template<int DG, int GP, typename T = FloatType>")
print("struct PSILagrangeImpl;")
for dg in [1, 3, 6, 8]:
    for dtype in ["float", "double"]:
        basisfunctions_in_lagrangepoints_cell(dg, 2, dtype)
    for dtype in ["float", "double"]:
        basisfunctions_in_lagrangepoints_cell(dg, 3, dtype)
print("template<int DG, int GP, typename T = FloatType>")
print(
    "const Eigen::Matrix<T, DG, GP*GP, Eigen::RowMajor> PSILagrange = PSILagrangeImpl<DG, GP, T>::value;"
)
print("")


print("\n\n//------------------------------ Basis Functions in Gauss Points (cell)\n")

print("template<int DG, int GP, typename T = FloatType>")
print("struct PSIImpl;")
for dg in [1, 3, 6, 8]:
    for dtype in ["float", "double"]:
        basisfunctions_in_gausspoints_cell(dg, 1, dtype)
    for dtype in ["float", "double"]:
        basisfunctions_in_gausspoints_cell(dg, 2, dtype)
    for dtype in ["float", "double"]:
        basisfunctions_in_gausspoints_cell(dg, 3, dtype)
    for dtype in ["float", "double"]:
        basisfunctions_in_gausspoints_cell(dg, 4, dtype)

print("template<int DG, int GP, typename T = FloatType>")
print(
    "const Eigen::Matrix<T, DG, GP*GP, Eigen::RowMajor> PSI = PSIImpl<DG, GP, T>::value;"
)
print("")


print("template<int DG, int GP, typename T = FloatType>")
print("struct PSIxImpl;")
print("template<int DG, int GP, typename T = FloatType>")
print("struct PSIyImpl;")
for dg in [1, 3, 6, 8]:
    for dtype in ["float", "double"]:
        basisfunctions_in_gausspoints_cell_gradient(dg, 1, dtype)
    for dtype in ["float", "double"]:
        basisfunctions_in_gausspoints_cell_gradient(dg, 2, dtype)
    for dtype in ["float", "double"]:
        basisfunctions_in_gausspoints_cell_gradient(dg, 3, dtype)
    for dtype in ["float", "double"]:
        basisfunctions_in_gausspoints_cell_gradient(dg, 4, dtype)
print("template<int DG, int GP, typename T = FloatType>")
print(
    "const Eigen::Matrix<T, DG, GP*GP, Eigen::RowMajor> PSIx = PSIxImpl<DG, GP, T>::value;"
)
print("template<int DG, int GP, typename T = FloatType>")
print(
    "const Eigen::Matrix<T, DG, GP*GP, Eigen::RowMajor> PSIy = PSIyImpl<DG, GP, T>::value;"
)
print("")


for dg in [1, 3, 6, 8]:
    for dtype in ["float", "double"]:
        integration_basisfunctions_in_gausspoints_cell(dg, 2, dtype)
    for dtype in ["float", "double"]:
        integration_basisfunctions_in_gausspoints_cell(dg, 3, dtype)


print("\n\n//------------------------------ Edge Basis Functions in Gauss Points\n")
print("template<int DG, int GP, typename T = FloatType>")
print("struct PSIeImpl;")
for dg in [1, 2, 3]:
    for dtype in ["float", "double"]:
        edge_basisfunctions_in_gausspoints(dg, dg, dtype)
    print("")
for dtype in ["float", "double"]:
    edge_basisfunctions_in_gausspoints(2, 3, dtype)

print("template<int DG, int GP, typename T = FloatType>")
print(
    "const Eigen::Matrix<T, DG, GP, Eigen::RowMajor> PSIe = PSIeImpl<DG, GP, T>::value;"
)
print("")

# Some output
print("\n} /* namespace Nextsim */\n")
print(f"#endif /* {include_guard} */")
