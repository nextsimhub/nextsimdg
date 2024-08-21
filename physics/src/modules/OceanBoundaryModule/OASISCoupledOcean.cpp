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

void OASISCoupledOcean::setMetadata(const ModelMetadata& metadata)
{
    // The parent class knows how to set the communicators and partitions
    OASISCoupled::setMetadata(metadata);

#ifdef USE_OASIS
    // OASIS defining variable

    /* Populate cplIdIn with return values from def_var, based on input from cplStringsIn. We assume
     * that the coupledIdIn enum starts at zero, is continuously numbered, and is in the same order
     * as cplStringsIn. Same story vor the input variables below. */
    cplIdIn.resize(cplStringsIn.size());
    for (int i = 0; i < cplStringsIn.size(); ++i) {
        OASIS_CHECK_ERR(oasis_c_def_var(
            &cplIdIn[i], cplStringsIn[i].c_str(), partitionID, bundleSize, OASIS_IN, OASIS_DOUBLE));
    }

    cplIdIn.resize(cplStringsOut.size());
    for (int i = 0; i < cplStringsOut.size(); ++i) {
        OASIS_CHECK_ERR(oasis_c_def_var(&cplIdOut[i], cplStringsOut[i].c_str(), partitionID,
            bundleSize, OASIS_IN, OASIS_DOUBLE));
    }

    // OASIS finalising definition
    OASIS_CHECK_ERR(oasis_c_enddef());
#else
    std::string message = __func__ + std::string(": OASIS support not compiled in.\n");
    throw std::runtime_error(message);
#endif
}

void OASISCoupledOcean::updateBefore(const TimestepTime& tst)
{
    // Directly set the array values
#ifdef USE_OASIS

    int kinfo;
    const int dimension0 = ModelArray::dimensions(ModelArray::Type::H)[0];
    const int dimension1 = ModelArray::dimensions(ModelArray::Type::H)[1];

    OASIS_CHECK_ERR(oasis_c_get(cplIdIn[couplingIdIn::SST], OASISTime, dimension0, dimension1,
        bundleSize, OASIS_DOUBLE, OASIS_COL_MAJOR, &sst[0], &kinfo));

    OASIS_CHECK_ERR(oasis_c_get(cplIdIn[couplingIdIn::SSS], OASISTime, dimension0, dimension1,
        bundleSize, OASIS_DOUBLE, OASIS_COL_MAJOR, &sss[0], &kinfo));

    OASIS_CHECK_ERR(oasis_c_get(cplIdIn[couplingIdIn::UOCEAN], OASISTime, dimension0, dimension1,
        bundleSize, OASIS_DOUBLE, OASIS_COL_MAJOR, &u[0], &kinfo));

    OASIS_CHECK_ERR(oasis_c_get(cplIdIn[couplingIdIn::VOCEAN], OASISTime, dimension0, dimension1,
        bundleSize, OASIS_DOUBLE, OASIS_COL_MAJOR, &v[0], &kinfo));

    // TODO: Implement ssh reading and passing to dynamics!
    //    OASIS_CHECK_ERR(oasis_c_get(cplIdIn[couplingIdIn::SSH], OASISTime, dimension0, dimension1,
    //                                bundleSize, OASIS_DOUBLE, OASIS_COL_MAJOR, &ssh[0], &kinfo));

    // TODO: Handle mld being an optional field
    mld = 10.;

    cpml = Water::cp * Water::rho * mld;

    overElements(
        std::bind(&OASISCoupledOcean::updateTf, this, std::placeholders::_1, std::placeholders::_2),
        TimestepTime());

    Module::getImplementation<IIceOceanHeatFlux>().update(tst);

#else
    std::string message = __func__ + std::string(": OASIS support not compiled in.\n");
    throw std::runtime_error(message);
#endif
}

void OASISCoupledOcean::updateAfter(const TimestepTime& tst)
{
#ifdef USE_OASIS
    // TODO: Add OASIS send-calls here
    int kinfo;
    const int dimension0 = ModelArray::dimensions(ModelArray::Type::H)[0];
    const int dimension1 = ModelArray::dimensions(ModelArray::Type::H)[1];

    // We still need the the outputs
    HField dummy;
    dummy.resize();
    dummy.setData(0.);
    OASIS_CHECK_ERR(oasis_c_put(cplIdOut[couplingIdOut::TAUX], OASISTime, dimension0, dimension1, 1,
        OASIS_DOUBLE, OASIS_COL_MAJOR, &dummy[0], OASIS_No_Restart, &kinfo));

    OASIS_CHECK_ERR(oasis_c_put(cplIdOut[couplingIdOut::TAUY], OASISTime, dimension0, dimension1, 1,
        OASIS_DOUBLE, OASIS_COL_MAJOR, &dummy[0], OASIS_No_Restart, &kinfo));

    OASIS_CHECK_ERR(oasis_c_put(cplIdOut[couplingIdOut::EMP], OASISTime, dimension0, dimension1, 1,
        OASIS_DOUBLE, OASIS_COL_MAJOR, &dummy[0], OASIS_No_Restart, &kinfo));

    OASIS_CHECK_ERR(oasis_c_put(cplIdOut[couplingIdOut::QSW], OASISTime, dimension0, dimension1, 1,
        OASIS_DOUBLE, OASIS_COL_MAJOR, &dummy[0], OASIS_No_Restart, &kinfo));

    OASIS_CHECK_ERR(oasis_c_put(cplIdOut[couplingIdOut::QNOSUN], OASISTime, dimension0, dimension1,
        1, OASIS_DOUBLE, OASIS_COL_MAJOR, &dummy[0], OASIS_No_Restart, &kinfo));

    OASIS_CHECK_ERR(oasis_c_put(cplIdOut[couplingIdOut::SFLX], OASISTime, dimension0, dimension1, 1,
        OASIS_DOUBLE, OASIS_COL_MAJOR, &dummy[0], OASIS_No_Restart, &kinfo));

    OASIS_CHECK_ERR(oasis_c_put(cplIdOut[couplingIdOut::TAUMOD], OASISTime, dimension0, dimension1,
        1, OASIS_DOUBLE, OASIS_COL_MAJOR, &dummy[0], OASIS_No_Restart, &kinfo));

    OASIS_CHECK_ERR(oasis_c_put(cplIdOut[couplingIdOut::CICE], OASISTime, dimension0, dimension1, 1,
        OASIS_DOUBLE, OASIS_COL_MAJOR, &dummy[0], OASIS_No_Restart, &kinfo));

    // Increment the "OASIS" time by the number of seconds in the time step
    updateOASISTime(tst);
#else
    std::string message = __func__ + std::string(": OASIS support not compiled in.\n");
    throw std::runtime_error(message);
#endif
}

} /* namespace Nextsim */
