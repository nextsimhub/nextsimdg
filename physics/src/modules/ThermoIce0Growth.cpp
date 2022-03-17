/*!
 * @file ThermoIce0Growth.cpp
 *
 * @date Mar 17, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/ThermoIce0Growth.hpp"
#include "include/ModelArray.hpp"

namespace Nextsim {

void ThermoIce0Growth::update(const TimePoint& tsInitialTime)
{
    // Loop over the entire HField domain
    for (size_t i = 0; i < ModelArray::size(ModelArray::HField); ++i) {
    }

}
} /* namespace Nextsim */
