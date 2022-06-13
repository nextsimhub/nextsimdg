/*!
 * @file DevGrid.cpp
 *
 * @date Dec 20, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/DevGrid.hpp"

#include <cstddef>
#include <vector>

namespace Nextsim {

const std::string DevGrid::structureName = "devgrid";
const std::string DevGrid::xDimName = "x";
const std::string DevGrid::yDimName = "y";
const std::string DevGrid::nIceLayersName = "nLayers";
const int DevGrid::nx = 10;

} /* namespace Nextsim */
