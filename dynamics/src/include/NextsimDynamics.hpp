/*!
 * @author  Thomas Richter <thomas.richter@ovgu.de>
 */

#ifndef __NEXTSIMDYNAMICS_HPP
#define __NEXTSIMDYNAMICS_HPP

#include "include/FloatType.hpp"

#include <array>
/*!
 * This file include all the definitions that are shared by nearly all nextsim-dynamics
 * files.
 */

namespace Nextsim {

constexpr FloatType pi = 3.14159265358979323846;

enum TENSOR_ELEMENT { I11, I12, I22, N_TENSOR_ELEMENTS };
const std::array<TENSOR_ELEMENT, N_TENSOR_ELEMENTS> tensorElements = { I11, I12, I22 };
typedef std::array<FloatType, N_TENSOR_ELEMENTS> SymmetricTensor;
/*!
 *
 * The different mappings used
 *
 * CARTESIAN: std. 2d cartesian map. No metric terms
 * SPHERICAL: mapping to spherical system in x=lon / y=lat
 *            ( R cos(lat)cos(lon), R cos(lat)sin(lon), R sin(lat) )
 */
enum COORDINATES { SPHERICAL, CARTESIAN };

/*!
 * Radius of the earth in [m]
 */
constexpr FloatType EarthRadius = 6371000.0;

template <typename T> constexpr T SQR(T x) { return x * x; }

/*!
 * Computes the correct number of stress-unknowns
 * depending on the CG Discretization
 *
 * CG1 = 3, CG2 = 8
 */
#define CG2DGSTRESS(CG) ((CG == 1 ? 3 : (CG == 2 ? 8 : -1)))

// Number of local CG dofs per element
#define CGDOFS(CG) ((CG == 1) ? 4 : 9)

/*!
 * number of Gauss points depending on DG unknowns
 * 1 for dG(0)   DG=1
 * 2 for dG(1)   DG=3 and dQ(1) DG=4
 * 3 for dG(2)   DG=6, DG=8 and dQ(2) DG=9
 */
#define GAUSSPOINTS(Q)                                                                             \
    (((Q == 9) || (Q == 8) || (Q == 6)) ? 9 : ((Q == 3) || (Q == 4)) ? 4 : (Q == 1) ? 1 : -1)
#define GAUSSPOINTS1D(Q)                                                                           \
    (((Q == 9) || (Q == 8) || (Q == 6)) ? 3 : ((Q == 3) || (Q == 4)) ? 2 : (Q == 1) ? 1 : -1)

/*!
 * Polynomial degree of the DG basis, per coordinate direction (x or y taken
 * separately). This is not the same as DG itself, which counts the number of
 * basis functions:
 *   DG=1 -> {1}                                   degree 0
 *   DG=3 -> {1, x, y}                             degree 1
 *   DG=4 -> {1, x, y, xy}                         degree 1
 *   DG=6 -> {1, x, y, xy, x^2, y^2}               degree 2
 *   DG=8 -> {1, x, y, xy, x^2, y^2, x^2y, xy^2}   degree 2
 *   DG=9 -> DG=8 plus x^2y^2                      degree 2
 */
#define DGDEGREE(DG)                                                                               \
    (((DG == 9) || (DG == 8) || (DG == 6)) ? 2 : ((DG == 3) || (DG == 4)) ? 1 : (DG == 1) ? 0 : -1)

/*!
 * Number of Gauss points per direction needed to integrate a polynomial of
 * degree D exactly. A rule with N points is exact up to degree 2N-1.
 */
#define NGAUSSFORDEGREE1D(D) (((D) + 2) / 2)

#define EDGEDOFS(DG) ((DG == 1) ? 1 : (((DG == 3) || (DG == 4)) ? 2 : 3))

}

#endif // #define  __NEXTSIMDYNAMICS
