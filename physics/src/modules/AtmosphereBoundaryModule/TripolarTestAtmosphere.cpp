/*!
 * @author  Einar Ólason <einar.olason@nersc.no>
 */

#include "include/TripolarTestAtmosphere.hpp"

#include "include/BenchmarkCoordinates.hpp"

namespace Nextsim {

TripolarTestAtmosphere::TripolarTestAtmosphere() = default;

void TripolarTestAtmosphere::setData(const ModelState::DataMap& ms)
{
    // Directly set the array values
    qiaAccessor.getHostRW() = 0.;
    dqia_dtAccessor.getHostRW() = 0.;
    qowAccessor.getHostRW() = 0.;
    sublAccessor.getHostRW() = 0.;
    snowAccessor.getHostRW() = 0.;
    rainAccessor.getHostRW() = 0.;
    evapAccessor.getHostRW() = 0.;
    penSWAccessor.getHostRW() = 0.;

    uwindAccessor.getHostRW() = 0.;
    HField& vwind = vwindAccessor.getHostRW();
    // Perform the rest of the calculation per element
    for (size_t j = 0; j < BenchmarkCoordinates::ny(); ++j) {
        for (size_t i = 0; i < BenchmarkCoordinates::nx(); ++i) {
            if (BenchmarkCoordinates::fx()(i, j) < 0.5)
                vwind(i, j) = 20;
            else
                vwind(i, j) = -20;
        }
    }
}

void TripolarTestAtmosphere::update(const TimestepTime& tst) { IAtmosphereBoundary::update(tst); }

} /* namespace Nextsim */
