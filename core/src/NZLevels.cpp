/*!
 * @file NZLevels.cpp
 *
 * @date 26 Jan 2023
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/NZLevels.hpp"

namespace Nextsim {

size_t NZLevels::nZLevels;

void NZLevels::set(size_t n) { nZLevels = n; }

size_t NZLevels::get() { return nZLevels; }

} /* namespace Nextsim */
