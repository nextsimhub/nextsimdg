/*!
 * @author  Tim Spain <timothy.spain@nersc.no>
 */

#include "include/ERA5Atmosphere.hpp"

#include "include/Finalizer.hpp"
#include "include/ModelMetadata.hpp"
#include "include/NetCDFForcings.hpp"
#include "include/NextsimModule.hpp"
#include "include/ParaGridIO.hpp"
#include "include/constants.hpp"

#include <ctime>
#include <string>
#include <tuple>
#include <vector>

#include <ncDim.h>
#include <ncFile.h>
#include <ncVar.h>

namespace Nextsim {

static const std::string pfx = "ERA5Atmosphere";
static const std::string dirKey = pfx + ".directory";
static const std::string checkOverrideKey = pfx + ".ignore_missing";

static const std::map<int, std::string> keyMap = {
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
        { keyMap.at(DIRPATH_KEY), fileDirectory() },
        { keyMap.at(CHECKOVERRIDE_KEY), checkOverride() },
    };
}

void ERA5Atmosphere::update(const TimestepTime& tst)
{
    // TODO: Get more authoritative names for the forcings
    std::vector<std::pair<std::string, ModelArray*>> forcings = {
            { "tair", &tair },
            { "dew2m", &tdew },
            { "pair", &pair },
            { "sw_in", &sw_in },
            { "lw_in", &lw_in },
            { "wind_speed", &wind },
            { "u", &uwind },
            { "v", &vwind },
    };

    for (std::pair<std::string, ModelArray*>& entry : forcings) {
        try {
            *entry.second = getData(entry.first, tst.start);
        } catch (const std::exception& e) {
            throw std::runtime_error(std::string(e.what()) + "\nFailed to read ERA5 forcing " + entry.first);
        }
    }
    snow = 0; // FIXME get snow data
    rain = 0; // FIXME get rain data

    // Convert temperatures from ERA5 (K) to nextSIM (˚C)
    tair += celsius(0.);
    tdew += celsius(0.);
    fluxImpl->update(tst);

    try {
        checkFields();
    } catch (const std::exception& e) {
        throw std::runtime_error("ERA5Atmosphere:update: " + std::string(e.what()));
    }
}

void ERA5Atmosphere::setDirectory(const std::string& dir) { fileDirectory() = dir; }
const std::string& ERA5Atmosphere::getDirectory() { return fileDirectory(); }
const std::string ERA5Atmosphere::addDirectory(const std::string& file)
{
    return getDirectory() + "/" + file;
}
void ERA5Atmosphere::setData(const ModelState::DataMap& ms)
{
    IAtmosphereBoundary::setData(ms);
    fluxImpl->setData(ms);
}

std::string e5FilenameFromYear(const std::string& era5Name, size_t year)
{
    std::string filename = "ERA5_" + era5Name + "_y" + std::to_string(year) + ".nc";
    return ERA5Atmosphere::addDirectory(filename);
}

const std::string era5FromNSName(const std::string& nsName)
{
    static const std::map<std::string, std::string> era5FromNS = {
            {"tair", "t2m"},
            {"dew2m", "d2m"},
            {"pair", "msl"},
            {"sw_in", "msdwswrf"},
            {"lw_in", "msdwlwrf"},
            {"u", "u10"},
            {"v", "v10"},
    };
    return era5FromNS.at(nsName);
}

using era5Buffer = Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;

std::pair<double, double> getLatitudeCoeffs(const std::string& filename)
{
    netCDF::NcFile ncFile(filename, netCDF::NcFile::read, netCDF::NcFile::nc4);
    netCDF::NcVar latitudeVar(ncFile.getVar("latitude"));
    size_t nLat = latitudeVar.getDim(0).getSize();
    era5Buffer latitudeBuff(nLat, 1);
    latitudeVar.getVar({0,}, {nLat,}, latitudeBuff.data());
    ncFile.close();
    return {latitudeBuff(1, 0) - latitudeBuff(0, 0), latitudeBuff(0, 0)};
}

std::pair<double, double> latitudeCoeffs(const std::string& key)
{
    static std::map<std::string, std::pair<double, double>> coeffMap;
    if (!coeffMap.count(key)) {
        coeffMap[key] = getLatitudeCoeffs(key);
    }
    return coeffMap.at(key);
}

ModelArray getERA5VarIndexData(const std::string& era5Name, size_t year, size_t tIndex, const ModelArray& modelLon, const ModelArray& modelLat)
{
    std::string filePath = e5FilenameFromYear(era5Name, year);
    era5Buffer e5data = NetCDFForcings::getFileIndexData(filePath, tIndex);

    // x and y positions in the source grid for each point on the target grid
    auto [dlat, lat0] = latitudeCoeffs(filePath);
    double dlon = std::fabs(dlat);
    ModelArray iFrac = modelLon / dlon;
    ModelArray jFrac = (modelLat - lat0) / dlat;
    return NetCDFForcings::maFromBuffer(e5data, iFrac, jFrac);
}

// time index from a TimePoint
size_t timeIndex(const TimePoint& time)
{
    return (time.doy() - 1) * 24 + time.hour();
}

// Time interpolation happens here
ModelArray getERA5VarTimeData(const std::string& era5Name, const TimePoint& time, const ModelArray& modelLon, const ModelArray& modelLat)
{
    ModelArray v1 = getERA5VarIndexData(era5Name, time.year(), timeIndex(time), modelLon, modelLat);
    TimePoint t2 = time + Duration(3600);
    ModelArray v2 = getERA5VarIndexData(era5Name, t2.year(), timeIndex(t2), modelLon, modelLat);
    double f = t2.minute() / 60.;
    return v2 * f + v1 * (1-f);
}

ModelArray ERA5Atmosphere::maHypot(const ModelArray& x, const ModelArray& y) const
{
    if (x.trueSize() != y.trueSize()) {
        throw std::runtime_error("maHypot: Cannot operate on differently sized ModelArrays.");
    }
    ModelArray wind = 0*x;
    ModelArray* pWind = &wind;
    ModelComponent::overElements([x, y, pWind](size_t i, const TimestepTime& tst) { (*pWind)[i] = std::sqrt(x[i]*x[i] + y[i]*y[i]);}, TimestepTime());
    return wind;
}

const ModelArray ERA5Atmosphere::getData(const std::string& nsName, const TimePoint& time) const
{
    const ModelMetadata& meta = ModelMetadata::getInstance();
    const ModelArray& modelLon = meta.longitude();
    const ModelArray& modelLat = meta.latitude();
    if (nsName == "wind_speed") {
        ModelArray u, v;
        u = getERA5VarTimeData("u10", time, modelLon, modelLat);
        v = getERA5VarTimeData("v10", time, modelLon, modelLat);
        return maHypot(u, v);
    } else {
        return getERA5VarTimeData(era5FromNSName(nsName), time, modelLon, modelLat);
    }
}

} /* namespace Nextsim */
