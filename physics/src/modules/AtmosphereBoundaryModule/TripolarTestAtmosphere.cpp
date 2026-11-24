/*!
 * @author  Einar Ólason <einar.olason@nersc.no>
 */

#include "include/TripolarTestAtmosphere.hpp"

#include "include/BenchmarkCoordinates.hpp"
#include "include/ModelArray.hpp"

namespace Nextsim {

TripolarTestAtmosphere::TripolarTestAtmosphere() = default;

void TripolarTestAtmosphere::setData(const ModelState::DataMap& ms)
{
    IAtmosphereBoundary::setData(ms);
    BenchmarkCoordinates::setData();

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

    // Use global indices so the split is based on the position in the
    // global domain, not in each rank's local subdomain
    const size_t i0 = ModelArray::definedDimensions.at(ModelArray::Dimension::X).start;
    const size_t nxGlobal = ModelArray::definedDimensions.at(ModelArray::Dimension::X).globalLength;
    // Perform the rest of the calculation per element
    for (size_t j = 0; j < BenchmarkCoordinates::ny(); ++j) {
        for (size_t i = 0; i < BenchmarkCoordinates::nx(); ++i) {
            if (static_cast<FloatType>(i0 + i) / nxGlobal < 0.5)
                vwind(i, j) = -20;
            else
                vwind(i, j) = -20;
        }
    }
}

void TripolarTestAtmosphere::update(const TimestepTime& tst) { IAtmosphereBoundary::update(tst); }

} /* namespace Nextsim */
