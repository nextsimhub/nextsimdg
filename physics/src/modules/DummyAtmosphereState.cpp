/*!
 * @file DummyAtmosphereState.cpp
 *
 * @date May 10, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/DummyAtmosphereState.hpp"

namespace Nextsim {
void DummyAtmosphereState::setData(const ModelState::DataMap&)
{
    tair = -1;
    tdew = -0.5;
    pair = 1e5;
    rmix = -1; // Use the dew point, rather than the mixing ratio
    sw_in = 0; // night
    lw_in = 311;
    snowfall = 0;
    evap_minus_precip = 0;

    windSpeed = 0;
}
} /* namespace Nextsim */
