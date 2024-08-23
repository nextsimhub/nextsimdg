/*!
 * @file BenchmarkOcean.cpp
 *
 * @date 19 Apr 2023
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/BenchmarkOcean.hpp"

#include "include/BenchmarkCoordinates.hpp"
#include "include/constants.hpp"

namespace Nextsim {

void BenchmarkOcean::setData(const ModelState::DataMap& ms)
{
    IOceanBoundary::setData(ms);
    BenchmarkCoordinates::setData();
    // The material parameters of the ocean are fixed to ensure no
    // thermodynamics occurs
    sst = -1.;
    sss = 32;
    tf = -1.8;
    mld = 10.;
    cpml = Water::rho * Water::cp * mld[0];
    qio = 0;
    ssh = 0.;

    // The time and length scales of the current generation function
    constexpr double L = 512000.; // Size of the domain in km
    constexpr double vMaxOcean = 0.01; // 1 cm/s in m/s
    constexpr double T = 8 * 86400; // 8 days in seconds

    // The currents are constant wrt time and space
    u = vMaxOcean * (2 * BenchmarkCoordinates::fy() - 1);
    v = vMaxOcean * (1 - 2 * BenchmarkCoordinates::fx());
}

} /* namespace Nextsim */
