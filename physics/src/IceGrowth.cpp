/*!
 * @file IceGrowth.cpp
 *
 * @date Mar 15, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/IceGrowth.hpp"

#include "include/constants.hpp"
#include "include/VerticalIceGrowthModule.hpp"

namespace Nextsim {

const std::map<int, std::string> keyMap = {
    { IceGrowth::VERTICAL_GROWTH_KEY, "VerticalIceModel" },
    { IceGrowth::LATERAL_GROWTH_KEY, "LateralIceModel" },
};

void IceGrowth::configure()
{
    // Configure the vertical and lateral growth modules
    iVertical = std::move(Module::getInstance<IVerticalIceGrowth>());
}

void IceGrowth::update(const TimestepTime& tsTime)
{
    iVertical->update(tsTime);

    // new ice formation
    overElements(std::bind(&IceGrowth::newIceFormation, this, std::placeholders::_1, std::placeholders::_2), tsTime);

}

void IceGrowth::newIceFormation(size_t i, const TimestepTime& tst)
{
    // Flux cooling the ocean from open water
    // TODO Add assimilation fluxes here
    double coolingFlux = qow[i];
    // Temperature change of the mixed layer during this timestep
    double deltaTml = -coolingFlux / mixedLayerBulkHeatCapacity[i] * tst.step;
    // Initial temperature
    double t0 = sst[i];
    // Freezing point temperature
    double tf0 = tf[i];
    // Final temperature
    double t1 = t0 + deltaTml;

    // deal with cooling below the freezing point
    if (t1 < tf0) {
        // Heat lost cooling the mixed layer to freezing point
        double sensibleFlux = (tf0 - t0) / deltaTml * coolingFlux;
        // Any heat beyond that is latent heat forming new ice
        double latentFlux = coolingFlux - sensibleFlux;

        qow[i] = sensibleFlux;
        newice[i]
            = latentFlux * tst.step * (1 - cice[i]) / (Ice::Lf * Ice::rho);
    }
}
} /* namespace Nextsim */
