/*
 * @file BasicIceOceanHeatFlux.cpp
 *
 * @date Oct 19, 2021
 * @author timpai
 */

#include "include/BasicIceOceanHeatFlux.hpp"

#include "include/ExternalData.hpp"
#include "include/PrognosticData.hpp"
#include "include/constants.hpp"

namespace Nextsim {

double flux(const PrognosticData& prog, const ExternalData& exter, PhysicsData& phys, NextsimPhysics& nsp)
{
    // The ice bottom temperature is the freezing point of the surface seawater
    //TODO: Implement the salinity dependent freezing temperature
    double iceBottomTemperature = -1.8;
    double tDiff = prog.seaSurfaceTemperature() - iceBottomTemperature;

    // Temperature difference multiplied by bulk heat capacity
    double bulkHeat = tDiff * Water::rho * Water::cp;

    // Transfer rate depends on the mixed layer depth an d a timescale. Here, it is the timestep
    return bulkHeat * exter.mixedLayerDepth() / prog.timestep();
}

} /* namespace Nextsim */
