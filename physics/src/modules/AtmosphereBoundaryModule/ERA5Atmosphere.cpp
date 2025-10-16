/*!
 * @author  Tim Spain <timothy.spain@nersc.no>
 */

#include "include/ERA5Atmosphere.hpp"

#include "include/Finalizer.hpp"
#include "include/ModelMetadata.hpp"
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
        } else if (entry.first == "time" || entry.first == "valid_time") {
        } else {
            dataVar = entry.second;
        }
    }
    std::vector<size_t> start = {tIndex, 0, 0};
    std::vector<size_t> count = {1, nLat, nLon};
    // Time dimension
    start[0] = tIndex;
    count[0] = 1;

    data.resize(nLon, nLat);

    dataVar.getVar(start, count, data.data());
    static const std::string offset_name = "add_offset";
    static const std::string scale_name = "scale_factor";
    auto dataAtts = dataVar.getAtts();
    double a = 1.;
    if (dataAtts.count(scale_name)){
        dataAtts.at(scale_name).getValues(&a);
    }
    double b = 0.;
    if (dataAtts.count(offset_name)) {
        dataAtts.at(offset_name).getValues(&b);
    }
    ncFile.close();

    data *= a;
    data += b;
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

    int nxe5 = buffer.rows();
    int nye5 = buffer.cols();
    double ptsPerDegree = nxe5 / 360;
    int nxsrc = nxe5 + 2;
    int nysrc = nye5 + 2;

    era5Buffer srcBuffer(nxsrc, nysrc);

    srcBuffer(Eigen::seq(1, Eigen::last-1), Eigen::seq(1, Eigen::last-1)) = buffer;
    // Wrap-around columns at the x edges
    srcBuffer(0, Eigen::seq(1, Eigen::last-1)) = buffer(Eigen::last, Eigen::all);
    srcBuffer(Eigen::last, Eigen::seq(1, Eigen::last-1)) = buffer(0, Eigen::all);
    // Duplicate rows at the y edges
    srcBuffer(Eigen::all, 0) = srcBuffer(Eigen::all, 1);
    srcBuffer(Eigen::all, Eigen::last) = srcBuffer(Eigen::all, Eigen::last - 1);

    // lambdas to translate latitude and longitude to (fractional) index in
    // srcBuffer, including wrap-around columns.
    auto xFromLon = [ptsPerDegree](double lon) { return ptsPerDegree * lon + 1; };
    auto yFromLat = [ptsPerDegree, nye5](double lat) { return ptsPerDegree * (90. - lat) + 1; };

    ModelArray maData(ModelArray::Type::H);
    maData.resize();
    for (size_t j = 0; j < nyma; ++j) {
        for (size_t i = 0; i < nxma; ++i) {
            double iFloat = xFromLon(destLon(i, j));
            double jFloat = yFromLat(destLat(i, j));
            int ilo = iFloat;
            int ihi = ilo + 1;
            int jlo = jFloat;
            int jhi = jlo + 1;
            double fx = 1 - (iFloat - ilo);
            double fy = 1 - (jFloat - jlo);
            maData(i, j) = fx * fy * srcBuffer(ilo, jlo) +
                    (1 - fx) * fy * srcBuffer(ihi, jlo) +
                    fx * (1 - fy) * srcBuffer(ilo, jhi) +
                    (1 - fx) * (1 - fy) * srcBuffer(ihi, jhi);
        }
    }
    return maData;
}

era5Buffer era5BufferHypot(const era5Buffer& x, const era5Buffer& y)
{
    return (x.square() + y.square()).sqrt();
}

const ModelArray ERA5Atmosphere::getData(const std::string& nsName, const TimePoint& time) const
{
    const ModelMetadata& meta = ModelMetadata::getInstance();
    const ModelArray& modelLon = meta.longitude();
    const ModelArray& modelLat = meta.latitude();
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
