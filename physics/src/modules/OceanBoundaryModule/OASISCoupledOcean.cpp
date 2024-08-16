/*!
 * @file OASISCoupledOcean.cpp
 *
 * @date Sep 26, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 * @author Einar Ólason <einar.olason@nersc.no>
 */

#include "include/OASISCoupledOcean.hpp"
#include "include/IIceOceanHeatFlux.hpp"
#include "include/Module.hpp"
#include "include/constants.hpp"

namespace Nextsim {
OASISCoupledOcean::OASISCoupledOcean()
    : IOceanBoundary()
{
}

void OASISCoupledOcean::setMetadata(const ModelMetadata& metadata)
{
    // The parent class knows how to set the communicators and partitions
    OASISCoupled::setMetadata(metadata);

#ifdef USE_OASIS
    // TODO: Insert OASIS def_var and end_def calls
#else
    std::string message = __func__ + std::string(": OASIS support not compiled in.\n");
    throw std::runtime_error(message);
#endif
}

void OASISCoupledOcean::updateBefore(const TimestepTime& tst)
{
    // Directly set the array values
#ifdef USE_OASIS
    // TODO: Replace this code with OASIS receive-calls
#else
    std::string message = __func__ + std::string(": OASIS support not compiled in.\n");
    throw std::runtime_error(message);
#endif
    sss = 32.;
    u = 0;
    v = 0;
    mld = 10.;
    double tf32 = -1.751; // Hand calculated from S = 32 using UNESCO
    tf = tf32;
    sst = tf32; // Tf == SST ensures that there is no ice-ocean heat flux
    cpml = Water::cp * Water::rho * mld;
    qio = 0.;

    Module::getImplementation<IIceOceanHeatFlux>().update(tst);
}

void OASISCoupledOcean::updateAfter(const TimestepTime& tst)
{
#ifdef USE_OASIS
    // TODO: Add OASIS send-calls here
#else
    std::string message = __func__ + std::string(": OASIS support not compiled in.\n");
    throw std::runtime_error(message);
#endif
}

} /* namespace Nextsim */
