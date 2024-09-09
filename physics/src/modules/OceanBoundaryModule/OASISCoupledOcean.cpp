/*!
 * @file OASISCoupledOcean.cpp
 *
 * @date 09 Sep 2024
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
#ifdef USE_OASIS
    // OASIS defining variable

    /* Populate the couplingId map with the id string and number pair. We need to do this seperately
     * for the input (get) and output (put) variables. */
    // TODO: coplingID should be a map of <std::string, std::pair> where the pair is idNumber and
    // pointer to the ModelArray
    for (std::string idString : cplStringsIn) {
        int idNumber;
        OASIS_CHECK_ERR(oasis_c_def_var(&idNumber, idString.c_str(), metadata.OASISPartitionId,
            bundleSize, OASIS_IN, OASIS_DOUBLE));
        couplingId[idString] = idNumber;
    }

    for (std::string idString : cplStringsOut) {
        int idNumber;
        OASIS_CHECK_ERR(oasis_c_def_var(&idNumber, idString.c_str(), metadata.OASISPartitionId,
            bundleSize, OASIS_OUT, OASIS_DOUBLE));
        couplingId[idString] = idNumber;
    }
#else
    throw std::runtime_error(std::string(__func__) + ": " + OASISError);
#endif
}

void OASISCoupledOcean::updateBefore(const TimestepTime& tst)
{
    // Directly set the array values
#ifdef USE_OASIS

    int kinfo;
    const int dimension0 = ModelArray::dimensions(ModelArray::Type::H)[0];
    const int dimension1 = ModelArray::dimensions(ModelArray::Type::H)[1];

    OASIS_CHECK_ERR(oasis_c_get(couplingId.at(SSTKey), OASISTime, dimension0, dimension1,
        bundleSize, OASIS_DOUBLE, OASIS_COL_MAJOR, &sst[0], &kinfo));

    OASIS_CHECK_ERR(oasis_c_get(couplingId.at(SSSKey), OASISTime, dimension0, dimension1,
        bundleSize, OASIS_DOUBLE, OASIS_COL_MAJOR, &sss[0], &kinfo));

    OASIS_CHECK_ERR(oasis_c_get(couplingId.at(UOceanKey), OASISTime, dimension0, dimension1,
        bundleSize, OASIS_DOUBLE, OASIS_COL_MAJOR, &u[0], &kinfo));

    OASIS_CHECK_ERR(oasis_c_get(couplingId.at(VOceanKey), OASISTime, dimension0, dimension1,
        bundleSize, OASIS_DOUBLE, OASIS_COL_MAJOR, &v[0], &kinfo));

    // TODO: Implement ssh reading and passing to dynamics!
    //    OASIS_CHECK_ERR(oasis_c_get(cplIdIn[couplingIdIn::SSHKey], OASISTime, dimension0,
    //    dimension1,
    //                                bundleSize, OASIS_DOUBLE, OASIS_COL_MAJOR, &ssh[0], &kinfo));

    if (couplingId.find(SSHKey) != couplingId.end()) {
        OASIS_CHECK_ERR(oasis_c_get(couplingId.at(VOceanKey), OASISTime, dimension0, dimension1,
            bundleSize, OASIS_DOUBLE, OASIS_COL_MAJOR, &mld[0], &kinfo));
    } else {
        mld = firstLayerDepth;
    }

    cpml = Water::cp * Water::rho * mld;

    overElements(
        std::bind(&OASISCoupledOcean::updateTf, this, std::placeholders::_1, std::placeholders::_2),
        TimestepTime());

    Module::getImplementation<IIceOceanHeatFlux>().update(tst);

#else
    throw std::runtime_error(std::string(__func__) + ": " + OASISError);
#endif
}

void OASISCoupledOcean::updateAfter(const TimestepTime& tst)
{
#ifdef USE_OASIS
    int kinfo;
    const int dimension0 = ModelArray::dimensions(ModelArray::Type::H)[0];
    const int dimension1 = ModelArray::dimensions(ModelArray::Type::H)[1];

    // TODO We still need the the actual data
    HField dummy;
    dummy.resize();
    dummy.setData(0.);
    OASIS_CHECK_ERR(oasis_c_put(couplingId.at(TauXKey), OASISTime, dimension0, dimension1, 1,
        OASIS_DOUBLE, OASIS_COL_MAJOR, &dummy[0], OASIS_No_Restart, &kinfo));

    OASIS_CHECK_ERR(oasis_c_put(couplingId.at(TauYKey), OASISTime, dimension0, dimension1, 1,
        OASIS_DOUBLE, OASIS_COL_MAJOR, &dummy[0], OASIS_No_Restart, &kinfo));

    OASIS_CHECK_ERR(oasis_c_put(couplingId.at(EMPKey), OASISTime, dimension0, dimension1, 1,
        OASIS_DOUBLE, OASIS_COL_MAJOR, &dummy[0], OASIS_No_Restart, &kinfo));

    OASIS_CHECK_ERR(oasis_c_put(couplingId.at(QSWKey), OASISTime, dimension0, dimension1, 1,
        OASIS_DOUBLE, OASIS_COL_MAJOR, &dummy[0], OASIS_No_Restart, &kinfo));

    OASIS_CHECK_ERR(oasis_c_put(couplingId.at(QNoSunKey), OASISTime, dimension0, dimension1, 1,
        OASIS_DOUBLE, OASIS_COL_MAJOR, &dummy[0], OASIS_No_Restart, &kinfo));

    OASIS_CHECK_ERR(oasis_c_put(couplingId.at(SFluxKey), OASISTime, dimension0, dimension1, 1,
        OASIS_DOUBLE, OASIS_COL_MAJOR, &dummy[0], OASIS_No_Restart, &kinfo));

    OASIS_CHECK_ERR(oasis_c_put(couplingId.at(TauModKey), OASISTime, dimension0, dimension1, 1,
        OASIS_DOUBLE, OASIS_COL_MAJOR, &dummy[0], OASIS_No_Restart, &kinfo));

    OASIS_CHECK_ERR(oasis_c_put(couplingId.at(CIceKey), OASISTime, dimension0, dimension1, 1,
        OASIS_DOUBLE, OASIS_COL_MAJOR, &dummy[0], OASIS_No_Restart, &kinfo));

    // Increment the "OASIS" time by the number of seconds in the time step
    updateOASISTime(tst);
#else
    throw std::runtime_error(std::string(__func__) + ": " + OASISError);
#endif
}

void OASISCoupledOcean::configure()
{
    SSTKey = Configured::getConfiguration(SSTConfigKey, SSTKeyDefault);
    SSSKey = Configured::getConfiguration(SSSConfigKey, SSSKeyDefault);
    UOceanKey = Configured::getConfiguration(UOceanConfigKey, UOceanKeyDefault);
    VOceanKey = Configured::getConfiguration(VOceanConfigKey, VOceanKeyDefault);
    SSHKey = Configured::getConfiguration(SSHConfigKey, SSHKeyDefault);
    MLDKey = Configured::getConfiguration(MLDConfigKey, MLDKeyDefault);
    TauXKey = Configured::getConfiguration(TauXConfigKey, TauXKeyDefault);
    TauYKey = Configured::getConfiguration(TauYConfigKey, TauYKeyDefault);
    TauModKey = Configured::getConfiguration(TauModConfigKey, TauModKeyDefault);
    EMPKey = Configured::getConfiguration(EMPConfigKey, EMPKeyDefault);
    QNoSunKey = Configured::getConfiguration(QNoSunConfigKey, QNoSunKeyDefault);
    QSWKey = Configured::getConfiguration(QSWConfigKey, QSWKeyDefault);
    SFluxKey = Configured::getConfiguration(SFluxConfigKey, SFluxKeyDefault);
    CIceKey = Configured::getConfiguration(CIceConfigKey, CIceKeyDefault);

    firstLayerDepth = Configured::getConfiguration(layerDepthConfigKey, FIRST_LAYER_DEPTH);
    if (Configured::getConfiguration(exchangeFirstLayerConfigKey, EXCHANGE_FIRST_LAYER)) {
        cplStringsIn.push_back(MLDKey);
    }
}

OASISCoupledOcean::HelpMap& OASISCoupledOcean::getHelpText(HelpMap& map, bool getAll)
{
    map[moduleName] = {
        { layerDepthConfigKey, ConfigType::NUMERIC, { "0", "∞" }, std::to_string(FIRST_LAYER_DEPTH),
            "m", "Depth of the first ocean model layer (if this is fixed)." },
        { exchangeFirstLayerConfigKey, ConfigType::BOOLEAN, { "true", "false" },
            std::to_string(EXCHANGE_FIRST_LAYER), "",
            "Use the thickness of the first ocean layer provided through the coupler" },
        { SSTConfigKey, ConfigType::STRING, {}, SSTKeyDefault, "",
            "The field name for sea surface temperature used in namcouple" },
        { SSSConfigKey, ConfigType::STRING, {}, SSSKeyDefault, "",
            "The field name for sea surface salinity used in namcouple" },
        { UOceanConfigKey, ConfigType::STRING, {}, UOceanKeyDefault, "",
            "The field name for ocean u-velocity used in namcouple" },
        { VOceanConfigKey, ConfigType::STRING, {}, VOceanKeyDefault, "",
            "The field name for ocean v-velocity used in namcouple" },
        { SSHConfigKey, ConfigType::STRING, {}, SSHKeyDefault, "",
            "The field name for sea surface height used in namcouple" },
        { MLDConfigKey, ConfigType::STRING, {}, MLDKeyDefault, "",
            "The field name for the thickness of the first ocean layer in namcouple (if that's "
            "defined)." },
        { TauXConfigKey, ConfigType::STRING, {}, TauXKeyDefault, "",
            "The field name for the x-component of ice-ocean stress in namcouple." },
        { TauYConfigKey, ConfigType::STRING, {}, TauYKeyDefault, "",
            "The field name for the y-component of ice-ocean stress in namcouple." },
        { TauModConfigKey, ConfigType::STRING, {}, TauModKeyDefault, "",
            "The field name for the modulus of ice-ocean stress in namcouple." },
        { EMPConfigKey, ConfigType::STRING, {}, EMPKeyDefault, "",
            "The field name for freshwater flux used in namcouple" },
        { QNoSunConfigKey, ConfigType::STRING, {}, QNoSunKeyDefault, "",
            "The field name for non-solar flux used in namcouple" },
        { QSWConfigKey, ConfigType::STRING, {}, QSWKeyDefault, "",
            "The field name for showrt-wave flux used in namcouple" },
        { QSWConfigKey, ConfigType::STRING, {}, QSWKeyDefault, "",
            "The field name for showrt-wave flux used in namcouple" },
        { SFluxConfigKey, ConfigType::STRING, {}, SFluxKeyDefault, "",
            "The field name for salt flux used in namcouple" },
        { CIceConfigKey, ConfigType::STRING, {}, CIceKeyDefault, "",
            "The field name for sea-ice concentration used in namcouple" },
    };
    return map;
}
OASISCoupledOcean::HelpMap& OASISCoupledOcean::getHelpRecursive(HelpMap& map, bool getAll)
{
    return getHelpText(map, getAll);
}

} /* namespace Nextsim */
