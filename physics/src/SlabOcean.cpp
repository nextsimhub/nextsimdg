/*!
 * @file SlabOcean.cpp
 *
 * @date 27 Jan 2023
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/SlabOcean.hpp"

#include "include/constants.hpp"

#include <map>
#include <string>

#include <oasis_c.h>

namespace Nextsim {

const std::string SlabOcean::sstSlabName = "sst_slab";
const std::string SlabOcean::sssSlabName = "sss_slab";
const double SlabOcean::defaultRelaxationTime = 30 * 24 * 60 * 60; // 30 days in seconds
int SlabOcean::date_cpl;
std::vector<int> SlabOcean::cplId;

// Configuration strings
static const std::string className = "SlabOcean";
static const std::string timeTName = "timeT";
static const std::string timeSName = "timeS";

template <>
const std::map<int, std::string> Configured<SlabOcean>::keyMap = {
    { SlabOcean::TIMET_KEY, className + "." + timeTName },
    { SlabOcean::TIMES_KEY, className + "." + timeSName },
};
void SlabOcean::configure()
{
    timeT = Configured::getConfiguration(keyMap.at(TIMET_KEY), defaultRelaxationTime);
    timeS = Configured::getConfiguration(keyMap.at(TIMES_KEY), defaultRelaxationTime);

    registerProtectedArray(ProtectedArray::SLAB_QDW, &qdw);
    registerProtectedArray(ProtectedArray::SLAB_FDW, &fdw);
    registerProtectedArray(ProtectedArray::SLAB_SST, &sstSlab);
    registerProtectedArray(ProtectedArray::SLAB_SSS, &sssSlab);
    registerProtectedArray(ProtectedArray::SLAB_CPL, &cplarr);
}

SlabOcean::HelpMap& SlabOcean::getHelpText(HelpMap& map, bool getAll)
{
    map[className] = {
        { keyMap.at(TIMET_KEY), ConfigType::NUMERIC, { "0", "∞" },
            std::to_string(defaultRelaxationTime), "s",
            "Relaxation time of the slab ocean to external temperature forcing." },
        { keyMap.at(TIMES_KEY), ConfigType::NUMERIC, { "0", "∞" },
            std::to_string(defaultRelaxationTime), "s",
            "Relaxation time of the slab ocean to external salinity forcing." },
    };
    return map;
};

void SlabOcean::setData(const ModelState::DataMap& ms)
{
    qdw.resize();
    fdw.resize();
    sstSlab.resize();
    sssSlab.resize();
    cplarr.resize();
    date_cpl = 0;

    // OASIS decomposition definition (global array on one MPI process)
    int part_params[OASIS_Serial_Params];
    int part_in;
    typedef ModelArray::Type Type;
    part_params[OASIS_Strategy] = OASIS_Serial;
    part_params[OASIS_Offset] = 0;
    part_params[OASIS_Length] = ModelArray::dimensions(Type::H)[0]*ModelArray::dimensions(Type::H)[1];
    OASIS_CHECK_ERR(oasis_c_def_partition(&part_in, OASIS_Serial_Params,
        				  part_params, part_params[OASIS_Length],
					  "part_in"));
    printf("OASIS grid exchange %d x %d \n ", ModelArray::dimensions(Type::H)[0],ModelArray::dimensions(Type::H)[1]);
    fflush(stdout);

    // OASIS defining variable
    int bundle_size;
    bundle_size = 1;
    cplId.resize(10);
    OASIS_CHECK_ERR(oasis_c_def_var(&cplId[1], "I_SSTSST", part_in, bundle_size,
  				    OASIS_IN, OASIS_DOUBLE));
    OASIS_CHECK_ERR(oasis_c_def_var(&cplId[2], "I_SSSal", part_in, bundle_size,
  				    OASIS_IN, OASIS_DOUBLE));
    OASIS_CHECK_ERR(oasis_c_def_var(&cplId[3], "I_OTaux1", part_in, bundle_size,
  				    OASIS_OUT, OASIS_DOUBLE));
    OASIS_CHECK_ERR(oasis_c_def_var(&cplId[4], "I_OTauy1", part_in, bundle_size,
  				    OASIS_OUT, OASIS_DOUBLE));
    OASIS_CHECK_ERR(oasis_c_def_var(&cplId[5], "IOEvaMPr", part_in, bundle_size,
  				    OASIS_OUT, OASIS_DOUBLE));
    OASIS_CHECK_ERR(oasis_c_def_var(&cplId[6], "I_QnsOce", part_in, bundle_size,
  				    OASIS_OUT, OASIS_DOUBLE));
    OASIS_CHECK_ERR(oasis_c_def_var(&cplId[7], "I_QsrOce", part_in, bundle_size,
  				    OASIS_OUT, OASIS_DOUBLE));
    OASIS_CHECK_ERR(oasis_c_def_var(&cplId[8], "I_SFLX", part_in, bundle_size,
  				    OASIS_OUT, OASIS_DOUBLE));
    OASIS_CHECK_ERR(oasis_c_def_var(&cplId[9], "I_TauMod", part_in, bundle_size,
  				    OASIS_OUT, OASIS_DOUBLE));
    OASIS_CHECK_ERR(oasis_c_def_var(&cplId[10], "IIceFrc", part_in, bundle_size,
  				    OASIS_OUT, OASIS_DOUBLE));
    // OASIS finalising definition
    OASIS_CHECK_ERR(oasis_c_enddef());
}

ModelState SlabOcean::getState() const
{
    return { {
                 { sstSlabName, sstSlab },
                 { sssSlabName, sssSlab },
             },
        {} };
}
ModelState SlabOcean::getState(const OutputLevel&) const { return getState(); }

std::unordered_set<std::string> SlabOcean::hFields() const { return { sstSlabName, sssSlabName }; }

void SlabOcean::update(const TimestepTime& tst)
{
    // OASIS sending field 
    int kinfo;
    typedef ModelArray::Type Type;

    cplarr = 0.;
    OASIS_CHECK_ERR(oasis_c_put(cplId[3], date_cpl, ModelArray::dimensions(Type::H)[0], ModelArray::dimensions(Type::H)[1],
                                1, OASIS_DOUBLE, OASIS_COL_MAJOR, &cplarr[0], OASIS_No_Restart, &kinfo));
    OASIS_CHECK_ERR(oasis_c_put(cplId[4], date_cpl, ModelArray::dimensions(Type::H)[0], ModelArray::dimensions(Type::H)[1],
                                1, OASIS_DOUBLE, OASIS_COL_MAJOR, &cplarr[0], OASIS_No_Restart, &kinfo));
    cplarr = emp;
    OASIS_CHECK_ERR(oasis_c_put(cplId[5], date_cpl, ModelArray::dimensions(Type::H)[0], ModelArray::dimensions(Type::H)[1],
                                1, OASIS_DOUBLE, OASIS_COL_MAJOR, &cplarr[0], OASIS_No_Restart, &kinfo));
    //WARNING  sw_in incorrect values
    cplarr = sw_in * 0.93 * ( 1 - cice );
    OASIS_CHECK_ERR(oasis_c_put(cplId[7], date_cpl, ModelArray::dimensions(Type::H)[0], ModelArray::dimensions(Type::H)[1],
                                1, OASIS_DOUBLE, OASIS_COL_MAJOR, &cplarr[0], OASIS_No_Restart, &kinfo));
    //cplarr = qow * ( cice - 1. ) + cplarr ;
    cplarr = qow * ( cice - 1 ) ;
    OASIS_CHECK_ERR(oasis_c_put(cplId[6], date_cpl, ModelArray::dimensions(Type::H)[0], ModelArray::dimensions(Type::H)[1],
                                1, OASIS_DOUBLE, OASIS_COL_MAJOR, &cplarr[0], OASIS_No_Restart, &kinfo));
    OASIS_CHECK_ERR(oasis_c_put(cplId[8], date_cpl, ModelArray::dimensions(Type::H)[0], ModelArray::dimensions(Type::H)[1],
                                1, OASIS_DOUBLE, OASIS_COL_MAJOR, &cplarr[0], OASIS_No_Restart, &kinfo));
    cplarr = v_air ;
    OASIS_CHECK_ERR(oasis_c_put(cplId[9], date_cpl, ModelArray::dimensions(Type::H)[0], ModelArray::dimensions(Type::H)[1],
                                1, OASIS_DOUBLE, OASIS_COL_MAJOR, &cplarr[0], OASIS_No_Restart, &kinfo));
    cplarr = cice;
    OASIS_CHECK_ERR(oasis_c_put(cplId[10], date_cpl, ModelArray::dimensions(Type::H)[0], ModelArray::dimensions(Type::H)[1],
                                1, OASIS_DOUBLE, OASIS_COL_MAJOR, &cplarr[0], OASIS_No_Restart, &kinfo));
    if (kinfo>0) {
    printf("sent field updated \n");
    }

    double dt = tst.step.seconds();
    // Slab SST update
    qdw = (sstExt - sst) * cpml / timeT;
    HField qioMean = qio * cice; // cice at start of TS, not updated
    HField qowMean = qow * (1 - cice); // 1- cice = open water fraction
    //sstSlab = sst - dt * (qioMean + qowMean - qdw) / cpml;
    // Slab SSS update
    HField arealDensity = cpml / Water::cp; // density times depth, or cpml divided by cp
    // This is simplified compared to the finiteelement.cpp calculation
    // Fdw = delS * mld * physical::rhow /(timeS*M_sss[i] - ddt*delS) where delS = sssSlab - sssExt
    fdw = (1 - sssExt / sss) * arealDensity / timeS;
    // ice volume change, both laterally and vertically
    HField deltaIceVol = newIce + deltaHice * cice;
    // change in snow volume due to melting (should be < 0)
    HField meltSnowVol = deltaSmelt * cice;
    // Mass per unit area after all the changes in water volume
    HField denominator
        = arealDensity - deltaIceVol * Ice::rho - meltSnowVol * Ice::rhoSnow - (emp - fdw) * dt;
    // Clamp the denominator to be at least 1 m deep, i.e. at least Water::rho kg m⁻²
    denominator.clampAbove(Water::rho);
    // Effective ice salinity is always less than or equal to the SSS
    HField effectiveIceSal = sss;
    effectiveIceSal.clampBelow(Ice::s);
    //sssSlab = sss
    //    + ((sss - effectiveIceSal) * Ice::rho * deltaIceVol // Change due to ice changes
    //          + sss * meltSnowVol
    //          + (emp - fdw) * dt) // snow melt, precipitation and nudging fluxes.
    //        / denominator;

    // OASIS receiving field 
    printf("OASIS ts %d \n ",date_cpl);
    fflush(stdout);
    OASIS_CHECK_ERR(oasis_c_get(cplId[1], date_cpl, ModelArray::dimensions(Type::H)[0],
                                ModelArray::dimensions(Type::H)[1], 1, OASIS_DOUBLE, OASIS_COL_MAJOR, &cplarr[0], &kinfo));
    if (kinfo>0) {
    sstSlab = cplarr;
    printf("SST updated \n");
    }
    printf("SSTOASIS %f %f\n ",sstSlab[100],sstSlab[1000]);
    fflush(stdout);
    OASIS_CHECK_ERR(oasis_c_get(cplId[2], date_cpl, ModelArray::dimensions(Type::H)[0],
                                ModelArray::dimensions(Type::H)[1], 1, OASIS_DOUBLE, OASIS_COL_MAJOR, &cplarr[0], &kinfo));
    if (kinfo>0) {
    sssSlab = cplarr;
    printf("SSS updated \n");
    }
    printf("SSSOASIS %f %f\n ",sssSlab[100],sssSlab[1000]);
    fflush(stdout);
    date_cpl+=tst.step.seconds();
}

} /* namespace Nextsim */
