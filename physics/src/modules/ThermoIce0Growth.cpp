/*!
 * @file ThermoIce0Growth.cpp
 *
 * @date Mar 17, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/ThermoIce0Growth.hpp"
#include "include/ModelArray.hpp"
#include "include/constants.hpp"

namespace Nextsim {

void ThermoIce0Growth::update(const TimePoint& tsInitialTime)
{
    // Loop over the entire HField domain
    for (size_t i = 0; i < ModelArray::size(ModelArray::Type::H); ++i) {
        calculateElement(i);
    }
}

void ThermoIce0Growth::calculateElement(size_t i)
{
    static const double bulkLHFusionSnow = Water::Lf * Ice::rhoSnow;

    double remainingFlux = *qic[i] - *qia[i];
    double snowMeltRate = std::min(-remainingFlux, 0.) / bulkLHFusionSnow;
    double snowSublRate = *sublim[i] / Ice::rhoSnow;
}
} /* namespace Nextsim */
