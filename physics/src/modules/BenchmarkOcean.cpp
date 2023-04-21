/*!
 * @file BenchmarkOcean.cpp
 *
 * @date 19 Apr 2023
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/BenchmarkOcean.hpp"

#include "include/constants.hpp"

namespace Nextsim {

void BenchmarkOcean::setData(const ModelState::DataMap&)
{
    // The material parameters of the ocean are fixed to ensure no
    // thermodynamics occurs
    sst = -1.;
    sss = 32;
    tf = -1.8;
    mld = 10.;
    cpml = Water::rho * Water::cp * mld[0];
    qio = 0;
}

void BenchmarkOcean::updateBefore(const TimestepTime& tst)
{
    // The time and length scales of the current generation function
    constexpr double L = 512000.; // Size of the domain in km
    constexpr double vMaxOcean = 0.01; // 1 cm/s in m/s
    constexpr double T = 8 * 86400; // 8 days in seconds


}
} /* namespace Nextsim */
