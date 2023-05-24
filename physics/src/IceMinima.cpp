/*!
 * @file IceMinima.cpp
 *
 * @date 8 Feb 2023
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/IceMinima.hpp"

namespace Nextsim {

const double IceMinima::hMinDefault = 0.01;
const double IceMinima::cMinDefault = 1e-12;

double IceMinima::hMin = hMinDefault;
double IceMinima::cMin = cMinDefault;

} /* namespace Nextsim */
