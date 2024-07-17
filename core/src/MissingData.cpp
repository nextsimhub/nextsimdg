/*!
 * @file MissingData.cpp
 *
 * @date Jul 17, 2024
 * @author Tim Spain <timothy.spain@nersc.no>
 * @author Einar Ólason <einar.olason@nersc.no>
 */

#include "include/MissingData.hpp"

namespace Nextsim {

const double MissingData::defaultValue = 1.7e38;
double MissingData::value = MissingData::defaultValue;

} /* namespace Nextsim */
