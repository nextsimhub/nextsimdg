/*!
 * @file MissingData.cpp
 *
 * @date Jun 14, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/MissingData.hpp"

namespace Nextsim {

const double MissingData::defaultValue = 1.7e38;
double MissingData::value = MissingData::defaultValue;

} /* namespace Nextsim */
