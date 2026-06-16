#ifndef __FLOATTYPE_HPP
#define __FLOATTYPE_HPP

namespace Nextsim {
/*!
 * The Float-Type used in NextSim for computations and outputs.
 * NEXTSIM_FLOAT_TYPE is defined in the main CMakeLists.txt based on the USE_SINGLE_PRECISION
 * option. Support for float is experimental and most useful when used with a Kokkos GPU backend.
 */
using FloatType = NEXTSIM_FLOAT_TYPE;

constexpr FloatType operator""_ft(long double v) { return FloatType(v); }

} /* namespace Nextsim */

#endif /* #define __FLOATTYPE_HPP */
