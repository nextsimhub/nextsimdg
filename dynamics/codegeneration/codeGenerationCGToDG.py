# generates c++ header files to define several static const matrices.
#
# here: projection between CG and DG vectors
#

import basisfunctions as bf
import gaussquadrature as gq
import numpy as np


# Compute the Projection matris realizing
#
# (dg, psi) = (cg, psi) in terms of dg = A * cg
#
# 'd' is the degree of the DG space
# dtype: C++ type name ("float" or "double")
def cg2dg_matrix(dg, cg, dtype):
    # print header
    if dg > 1:
        print(
            "static const Eigen::Matrix<{0}, {1}, {2}, Eigen::RowMajor> CG{3}_to_DG{4}_{5} =".format(
                dtype, dg, bf.cgdofs(cg), cg, dg, dtype
            )
        )
        print(
            "\t(Eigen::Matrix<{0}, {1}, {2}, Eigen::RowMajor>() <<".format(
                dtype, dg, bf.cgdofs(cg)
            )
        )
    else:
        print(
            "static const Eigen::Matrix<{0}, {1}, {2}> CG{3}_to_DG{4}_{5} =".format(
                dtype, dg, bf.cgdofs(cg), cg, dg, dtype
            )
        )
        print("\t(Eigen::Matrix<{0}, {1}, {2}>() <<".format(dtype, dg, bf.cgdofs(cg)))

    for dgi in range(dg):
        for cgi in range(bf.cgdofs(cg)):
            xxx = 0
            for gx in range(3):
                for gy in range(3):
                    X = gq.gausspoints[2][gx]
                    Y = gq.gausspoints[2][gy]
                    xxx = xxx + gq.gaussweights[2][gx] * gq.gaussweights[2][
                        gy
                    ] * bf.CGbasisfunction(cg, cgi, X, Y) * bf.dgbasis(dgi, X, Y)

            print(xxx * bf.inversemass[dgi], end="")
            if (cgi < bf.cgdofs(cg) - 1) or dgi < dg - 1:
                print(",", end="")
            else:
                print(").finished();")


# Compute the Projection matris realizing
#
# (dg, psi) = (d_X/Y cg, psi) in terms of dg = A_dX/Y * cg
#
# 'd' is the degree of the DG space
# dtype: C++ type name ("float" or "double")
def cg2dg_dxy_matrix(dg, cg, dXY, dtype):
    # print header
    if dg > 1:
        print(
            "static const Eigen::Matrix<{0}, {1}, {2}, Eigen::RowMajor> CG{3}_to_DG{4}_d{5}_{6} =".format(
                dtype, dg, bf.cgdofs(cg), cg, dg, dXY, dtype
            )
        )
        print(
            "\t(Eigen::Matrix<{0}, {1}, {2}, Eigen::RowMajor>() <<".format(
                dtype, dg, bf.cgdofs(cg)
            )
        )
    else:
        print(
            "static const Eigen::Matrix<{0}, {1}, {2}> CG{3}_to_DG{4}_d{5}_{6} =".format(
                dtype, dg, bf.cgdofs(cg), cg, dg, dXY, dtype
            )
        )
        print("\t(Eigen::Matrix<{0}, {1}, {2}>() <<".format(dtype, dg, bf.cgdofs(cg)))

    for dgi in range(dg):
        for cgi in range(bf.cgdofs(cg)):
            xxx = 0
            for gx in range(3):
                for gy in range(3):
                    X = gq.gausspoints[2][gx]
                    Y = gq.gausspoints[2][gy]
                    if dXY == "X":
                        xxx = xxx + gq.gaussweights[2][gx] * gq.gaussweights[2][
                            gy
                        ] * bf.CGbasisfunction_dX(cg, cgi, X, Y) * bf.dgbasis(dgi, X, Y)
                    elif dXY == "Y":
                        xxx = xxx + gq.gaussweights[2][gx] * gq.gaussweights[2][
                            gy
                        ] * bf.CGbasisfunction_dY(cg, cgi, X, Y) * bf.dgbasis(dgi, X, Y)
                    else:
                        print("Direction", dXY, "not known")
                        raise AssertionError

            print(xxx * bf.inversemass[dgi], end="")
            if (cgi < bf.cgdofs(cg) - 1) or dgi < dg - 1:
                print(",", end="")
            else:
                print(").finished();")


# Compute the Projection matris realizing
#
# (dg, psi) = (d_X/Y cg, psi) in terms of dg = A_dX/Y * cg
#
# 'd' is the degree of the DG space
# dtype: C++ type name ("float" or "double")
def dg_cg_dxy_matrix(dg, cg, dXY, dtype):
    # print header
    if dg > 1:
        print(
            "static const Eigen::Matrix<{0}, {1}, {2}, Eigen::RowMajor> DG{3}_CG{4}_d{5}_{6} =".format(
                dtype, bf.cgdofs(cg), dg, dg, cg, dXY, dtype
            )
        )
        print(
            "\t(Eigen::Matrix<{0}, {1}, {2}, Eigen::RowMajor>() <<".format(
                dtype, bf.cgdofs(cg), dg
            )
        )
    else:
        print(
            "static const Eigen::Matrix<{0}, {1}, {2}> DG{3}_CG{4}_d{5}_{6} =".format(
                dtype, bf.cgdofs(cg), dg, dg, cg, dXY, dtype
            )
        )
        print("\t(Eigen::Matrix<{0}, {1}, {2}>() <<".format(dtype, bf.cgdofs(cg), dg))

    for cgi in range(bf.cgdofs(cg)):
        for dgi in range(dg):
            xxx = 0
            for gx in range(3):
                for gy in range(3):
                    X = gq.gausspoints[2][gx]
                    Y = gq.gausspoints[2][gy]
                    if dXY == "X":
                        xxx = xxx + gq.gaussweights[2][gx] * gq.gaussweights[2][
                            gy
                        ] * bf.CGbasisfunction_dX(cg, cgi, X, Y) * bf.dgbasis(dgi, X, Y)
                    elif dXY == "Y":
                        xxx = xxx + gq.gaussweights[2][gx] * gq.gaussweights[2][
                            gy
                        ] * bf.CGbasisfunction_dY(cg, cgi, X, Y) * bf.dgbasis(dgi, X, Y)
                    else:
                        print("Direction", dXY, "not known")
                        raise AssertionError

            print(xxx * bf.inversecg[cg - 1][cgi], end="")
            if (cgi < bf.cgdofs(cg) - 1) or dgi < dg - 1:
                print(",", end="")
            else:
                print(").finished();")


### Main

# make sure that Guass quadrature is correct
gq.sanitycheck_gauss()

include_guard = "__CODEGENERATIONCGTODG_HPP"
# Some output
print(f"#ifndef {include_guard}")
print(f"#define {include_guard}")
print("\n")
print("// Automatically generated by codegeneration/project_cg_dg.py")
print("//")
print("// Generates the matrices CG2toDG[dg]")
print("// - Realizes the projection of a CG2 vector into the DG[dg] space")
print("//   dg = CG2toDG[dg] * cg")
print("")

print('#include "NextsimDynamics.hpp"')
print("namespace Nextsim {")

print("\n\n//------------------------------ CGtoDG\n")
for dg in [1, 3, 6, 8]:
    for cg in [1, 2]:
        for dtype in ["float", "double"]:
            cg2dg_matrix(dg, cg, dtype)
            print("")

print("\n")
print("// Generates the matrices CG2toDG[dg]_dX and CG2toDG[dg]_dY")
print(
    "// - Realizes the projection of the derivative of a CG2 vector into the DG[dg] space"
)
print("//   dg = CG2toDG[dg] * cg")
print("")


print("//------------------------------ CG2toDG\n")
for dg in [1, 3, 6, 8]:
    for cg in [1, 2]:
        for dtype in ["float", "double"]:
            cg2dg_dxy_matrix(dg, cg, "X", dtype)
            cg2dg_dxy_matrix(dg, cg, "Y", dtype)
            print("")

print("\n")
print("// Generates the matrices DG1_CG2_dX/Y for")
print("// adding the DG1 - Stress tensor to the CG2 equation, e.g. for (S, nabla Phi)")
print("// computed as += DG1_CG2_dX/Y * S11/12/22")
print("")


print("//------------------------------ CG2toDG\n")
for dg in [1, 3, 6, 8]:
    for cg in [1, 2]:
        for dtype in ["float", "double"]:
            dg_cg_dxy_matrix(dg, cg, "X", dtype)
            dg_cg_dxy_matrix(dg, cg, "Y", dtype)
            print("")


# Some output
print("\n} /* namespace Nextsim */\n")
print(f"#endif /* {include_guard} */")
