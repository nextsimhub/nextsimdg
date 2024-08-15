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
    // TODO: Insert OASIS initialisation calls here
    /*
     * These three may need to be global to the class
     * int compID;
     * MPI_Comm localComm;
     * MPI_Comm coupledComm;
     *
     * std::string compName = "nextsim";
     * std::string functionName = getName()+"::setMetadata";
     * std::string message = "couldn't initialise component";
     * std::string file = "OASISCoupledOcean.cpp";
     * const bool coupled = true;
     * if ( ! oasis_c_init_comp_with_comm(&compID, &compName.c_str(), coupled, metadata.mpiComm) )
     *     oasis_c_abort(compID, &functionName.c_str(), &message.c_str(), &file.c_str(), 36);
     *
     * if ( ! oasis_c_get_localcomm(&localComm) )
     *     oasis_c_abort(compID, &functionName.c_str(), &message.c_str(), &file.c_str(), 39);
     *
     * if ( ! oasis_c_create_couplcomm(metadata.mpiMyRank, &localComm, &coupledComm) )
     *     oasis_c_abort(compID, &functionName.c_str(), &message.c_str(), &file.c_str(), 42);
     *
     * ... and then def_partition, def_var and end_def calls
     *
     * But actually, everything except def_var (and the following end_def) can be seperated out in
     * its own OASISCoupled class, which OASISCoupledOcean and OASISCoupledAtmosphere inherit from
     * (and OASISCoupledWaves will also inherit from).
     */
}

void OASISCoupledOcean::updateBefore(const TimestepTime& tst)
{
    // Directly set the array values
    // TODO: Replace this code with OASIS receive-calls
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
    // TODO: Add OASIS send-calls here
}

OASISCoupledOcean::~OASISCoupledOcean()
{
        // TODO: Insert OASIS finalise call(s) here
}

} /* namespace Nextsim */
