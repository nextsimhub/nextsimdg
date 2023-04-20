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

}
} /* namespace Nextsim */
