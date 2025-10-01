/*!
 * @author  Tim Spain <timothy.spain@nersc.no>
 */

#include "include/ERA5Atmosphere.hpp"

#include "include/Finalizer.hpp"
#include "include/NextsimModule.hpp"
#include "include/ParaGridIO.hpp"

#include <ctime>

namespace Nextsim {

std::string ERA5Atmosphere::filePath;

static const std::string pfx = "ERA5Atmosphere";
static const std::string fileKey = pfx + ".file";
static const std::string dirKey = pfx + ".directory";
static const std::string checkOverrideKey = pfx + ".ignore_missing";

static const std::map<int, std::string> keyMap = {
    { ERA5Atmosphere::FILEPATH_KEY, fileKey },
    { ERA5Atmosphere::DIRPATH_KEY, dirKey },
    { ERA5Atmosphere::CHECKOVERRIDE_KEY, checkOverrideKey },
};

ERA5Atmosphere::ERA5Atmosphere()
    : fluxImpl(nullptr)
    , tair(ModelArray::Type::H, { -100, 100 })
    , tdew(ModelArray::Type::H, { -100, 100 })
    , pair(ModelArray::Type::H, { 500e2, 2000e2 })
    , sw_in(ModelArray::Type::H, { -1e-6, 1e4 })
    , lw_in(ModelArray::Type::H, { -1e-6, 1e4 })
    , wind(ModelArray::Type::H, { 0, 100 })
{
    getStore().registerArray(Protected::T_AIR, &tair, RO);
    getStore().registerArray(Protected::DEW_2M, &tdew, RO);
    getStore().registerArray(Protected::P_AIR, &pair, RO);
    getStore().registerArray(Protected::SW_IN, &sw_in, RO);
    getStore().registerArray(Protected::LW_IN, &lw_in, RO);
    getStore().registerArray(Protected::WIND_SPEED, &wind, RO);
}

ConfigurationHelp::HelpMap& ERA5Atmosphere::getHelpRecursive(HelpMap& map, bool getAll)
{
    map[pfx] = {
        { fileKey, ConfigType::STRING, {}, "", "",
            "Path to the processed NetCDF file providing the ERA5 forcings." },
        { dirKey, ConfigType::STRING, {}, "", "",
            "Path to the directory containing the ERA5 NetCDF files." },
        { dirKey, ConfigType::BOOLEAN, {"true", "false"}, "false", "",
            "Continue without error even if required ERA5 files are missing." },
    };
    Module::getHelpRecursive<IFluxCalculation>(map, getAll);

    return map;
}

void ERA5Atmosphere::configure()
{
    Finalizer::registerUnique(Module::finalize<IFluxCalculation>);

    filePath = Configured::getConfiguration(keyMap.at(FILEPATH_KEY), std::string());
    fileDirectory() = Configured::getConfiguration(keyMap.at(DIRPATH_KEY), fileDirectory());
    checkOverride() = Configured::getConfiguration(keyMap.at(CHECKOVERRIDE_KEY), false);

    fluxImpl = &Module::getImplementation<IFluxCalculation>();
    tryConfigure(fluxImpl);

    addChecks({
        { "tair", &tair },
        { "tdew", &tdew },
        { "pair", &pair },
        { "sw_in", &sw_in },
        { "lw_in", &lw_in },
        { "wind", &wind },
    });
}

ConfigMap ERA5Atmosphere::getConfiguration() const
{
    return {
        { keyMap.at(FILEPATH_KEY), filePath },
    };
}

void ERA5Atmosphere::update(const TimestepTime& tst)
{
    // TODO: Get more authoritative names for the forcings
    std::set<std::string> forcings
        = { "tair", "dew2m", "pair", "sw_in", "lw_in", "wind_speed", "u", "v" };

    ModelState state = ParaGridIO::readForcingTimeStatic(forcings, tst.start, filePath);
    tair = state.data.at("tair");
    tdew = state.data.at("dew2m");
    pair = state.data.at("pair");
    sw_in = state.data.at("sw_in");
    lw_in = state.data.at("lw_in");
    wind = state.data.at("wind_speed");
    uwind = state.data.at("u");
    vwind = state.data.at("v");
    snow = 0; // FIXME get snow data
    rain = 0; // FIXME get rain data

    fluxImpl->update(tst);

    try {
        checkFields();
    } catch (const std::exception& e) {
        throw std::runtime_error("ERA5Atmosphere:update: " + std::string(e.what()));
    }
}

void ERA5Atmosphere::setFilePath(const std::string& filePathIn) { filePath = filePathIn; }

void ERA5Atmosphere::setData(const ModelState::DataMap& ms)
{
    IAtmosphereBoundary::setData(ms);
    fluxImpl->setData(ms);
}

const std::string ERA5Atmosphere::filename(const std::string& era5Name, const TimePoint& time)
{
    std::tm* cTime = time.gmtime();
    int year = cTime->tm_year;
    return "ERA5_" + era5Name + "_y" + std::to_string(year) + ".nc";
}

const std::string ERA5Atmosphere::era5FromNSName(const std::string& nsName)
{
    static const std::map<std::string, std::string> era5FromNS = {
            {"tair", "t2m"},
            {"tdew", "d2m"},
            {"pair", "msl"},
            {"sw_in", "msdwswrf"},
            {"lw_in", "msdwlwrf"},
    };
    return era5FromNS.at(nsName);
}

const ModelArray ERA5Atmosphere::getData(const std::string& nsName, const TimePoint& time)
{
    if (nsName == "wind_speed") {
        era5Buffer u, v;
        u = dataBuffer(filename("u10", time), time);
        v = dataBuffer(filename("v10", time), time);
        return maFromERA5Buffer((u.square() + v.square()).sqrt());
    } else {
        return maFromERA5Buffer(dataBuffer(filename(era5FromNSName(nsName), time), time));
    }
}
} /* namespace Nextsim */
