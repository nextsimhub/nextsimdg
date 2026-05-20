/*!
 * @author Robert Jendersie <robert.jendersie@ovgu.de>
 */

#ifndef KERNEL_ALTERNATIVES_HPP
#define KERNEL_ALTERNATIVES_HPP

#include "FloatType.hpp"

#ifdef USE_KOKKOS
#include "../kokkos/include/KokkosUtils.hpp"
#include <Kokkos_Core.hpp>
#else
#include <algorithm>
#include <cmath>
#endif

namespace Nextsim {

// namespace Utils contains standard math functions (sin, sqrt, min, ...)
#ifdef USE_KOKKOS
namespace Utils = Kokkos;
#else
namespace Utils = std;
#endif

// KERNEL_IMPL_FUNCTION annotates functions that need to be callable inside a kernel
#ifdef USE_KOKKOS
#define KERNEL_IMPL_FUNCTION KOKKOS_FUNCTION
#else
// no special annotation needed
#define KERNEL_IMPL_FUNCTION
#endif

// KERNEL_LAMBDA capture list of a lambda used in a overElements call
#ifdef USE_KOKKOS
#define OVER_ELEMENTS_LAMBDA KOKKOS_LAMBDA
#define OVER_ELEMENTS_CLASS_LAMBDA KOKKOS_CLASS_LAMBDA
using ElementIndex = DeviceIndex;
#else
#define OVER_ELEMENTS_LAMBDA [&]
#define OVER_ELEMENTS_CLASS_LAMBDA [&]
using ElementIndex = size_t;
#endif

// Execution space that is used to signal async operations.
#ifdef USE_KOKKOS
using DefaultExecutionSpace = Kokkos::DefaultExecutionSpace;
#else
struct ExecutionSpaceDummy { };
using DefaultExecutionSpace = ExecutionSpaceDummy;
#endif

} /* namespace Nextsim */

#endif /* KERNEL_ALTERNATIVES_HPP */
