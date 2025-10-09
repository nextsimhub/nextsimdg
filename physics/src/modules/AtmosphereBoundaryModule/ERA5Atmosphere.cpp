/*!
 * @author  Tim Spain <timothy.spain@nersc.no>
 */

#include "include/ERA5Atmosphere.hpp"

#include "include/Finalizer.hpp"
#include "include/NextsimModule.hpp"
#include "include/ParaGridIO.hpp"

#include <ctime>

#include <ncDim.h>
#include <ncFile.h>
#include <ncVar.h>

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
//    return ERA5Atmosphere::addDirectory(filename);
    return ERA5Atmosphere::addDirectory(filename);
}

const std::string era5FromNSName(const std::string& nsName)
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

using era5Buffer = Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

era5Buffer getFileIndexData(const std::string& filename, size_t tIndex)
{
    era5Buffer data;
    netCDF::NcFile ncFile(filename, netCDF::NcFile::read, netCDF::NcFile::nc4);
    netCDF::NcVar dataVar;

    size_t nLat;
    size_t nLon;
    // There should be 4 vars, longitude, latitude, time and the data variable
    for (auto entry : ncFile.getVars()) {
        if (entry.first == "longitude") {
            nLon = entry.second.getDim(0).getSize();
        } else if (entry.first == "latitude") {
            nLat = entry.second.getDim(0).getSize();
        } else if (entry.first == "time") {
        } else {
            dataVar = entry.second;
        }
    }
    std::vector<size_t> start = {tIndex, 0, 0};
    std::vector<size_t> count = {1, nLat, nLon};
    // Time dimension
    start[0] = tIndex;
    count[0] = 1;

    data.resize(nLat, nLon);

    dataVar.getVar(start, count, data.data());

    ncFile.close();

    return data;
}

era5Buffer getVarIndexData(const std::string& era5Name, size_t year, size_t tIndex)
{
    return getFileIndexData(e5FilenameFromYear(era5Name, year), tIndex);
}

// time index from the C tm struct
size_t timeIndexFromTM(const std::tm* tm)
{
    return tm->tm_yday * 24 + tm->tm_hour;
}


// Time interpolation happens here
era5Buffer getVarTimeData(const std::string& era5Name, const TimePoint& time)
{
    std::tm* tm1 = time.gmtime();
    static size_t epochYear = 1900;

    era5Buffer v1 = getVarIndexData(era5Name, tm1->tm_year + epochYear, timeIndexFromTM(tm1));
    TimePoint t2 = time + Duration(3600);
    std::tm* tm2 = t2.gmtime();
    era5Buffer v2 = getVarIndexData(era5Name, tm2->tm_year + epochYear, timeIndexFromTM(tm1));
    double f = tm2->tm_min / 60.;
    return v2 * f + v1 * (1-f);
}

ModelArray maFromERA5Buffer(const era5Buffer& buffer, const HField& destLon, const HField& destLat)
{
    int nxma = destLon.dimensions()[0];
    int nyma = destLon.dimensions()[1];

    int nxe5 = 360*4;
    int nye5 = buffer.rows() / nxe5;

    int nxsrc = nxe5 + 2;
    int nysrc = nye5 + 2;

    era5Buffer srcBuffer(nxsrc, nysrc);

    return ModelArray(ModelArray::Type::H);
}

era5Buffer era5BufferHypot(const era5Buffer& x, const era5Buffer& y)
{
    return (x.square() + y.square()).sqrt();
}

const ModelArray ERA5Atmosphere::getData(const std::string& nsName, const TimePoint& time) const
{
    if (nsName == "wind_speed") {
        era5Buffer u, v;
        u = getVarTimeData("u10", time);
        v = getVarTimeData("v10", time);
        return maFromERA5Buffer(era5BufferHypot(u, v), modelLon, modelLat);
    } else {
        return maFromERA5Buffer(getVarTimeData(era5FromNSName(nsName), time), modelLon, modelLat);
    }
}

} /* namespace Nextsim */
