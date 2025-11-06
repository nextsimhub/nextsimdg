/*!
 * @author  Tim Spain <timothy.spain@nersc.no>
 */

#include "include/TOPAZOcean.hpp"

#include "include/Finalizer.hpp"
#include "include/NetCDFForcings.hpp"
#include "include/NextsimModule.hpp"
#include "include/ParaGridIO.hpp"
#include "include/constants.hpp"
#include "include/gridNames.hpp"

namespace Nextsim {

std::string TOPAZOcean::filePath;

static const std::string pfx = "TOPAZOcean";
static const std::string fileKey = pfx + ".file";

static const std::map<int, std::string> keyMap = {
    { TOPAZOcean::FILEPATH_KEY, fileKey },
};

TOPAZOcean::TOPAZOcean()
    : sstExt(ModelArray::Type::H, { -5, 50 })
    , sssExt(ModelArray::Type::H, { 0, 50 })
    , slabOcean(m_couplingArrays)
{
}

ConfigurationHelp::HelpMap& TOPAZOcean::getHelpRecursive(HelpMap& map, bool getAll)
{
    map[pfx] = {
        { fileKey, ConfigType::STRING, {}, "", "",
            "Path to the processed NetCDF file providing the TOPAZ forcings." },
    };

    return map;
}

void TOPAZOcean::configure()
{
    Finalizer::registerUnique(Module::finalize<IIceOceanHeatFlux>);
    Finalizer::registerUnique(Module::finalize<IFreezingPoint>);

    pIOHeatFlux = std::move(Module::getInstance<IIceOceanHeatFlux>());
    tryConfigure(*pIOHeatFlux);

    pFreezingPoint = std::move(Module::getInstance<IFreezingPoint>());
    tryConfigure(*pFreezingPoint);

    filePath = Configured::getConfiguration(keyMap.at(FILEPATH_KEY), std::string());

    slabOcean.configure();

    getStore().registerArray(Protected::EXT_SST, &sstExt, RO);
    getStore().registerArray(Protected::EXT_SSS, &sssExt, RO);
}

ConfigMap TOPAZOcean::getConfiguration() const { return { { keyMap.at(FILEPATH_KEY), filePath } }; }

void TOPAZOcean::updateBefore(const TimestepTime& tst)
{
    std::set<std::string> forcings = { sstName, sssName, mldName, uName, vName, sshName };

    ModelState state = ParaGridIO::readForcingTimeStatic(forcings, tst.start, filePath);
    sstExt = state.data.at(sstName);
    sssExt = state.data.at(sssName);
    mld = state.data.at(mldName);
    u = state.data.at(uName);
    v = state.data.at(vName);
    if (state.data.count(sshName)) {
        ssh = state.data.at(sshName);
    } else {
        ssh = 0.;
    }

    cpml = Water::rhoOcean * Water::cp * mld;
    overElements([this](size_t i, const TimestepTime& tsTime) { this->updateTf(i, tsTime); }, tst);

    Module::getImplementation<IIceOceanHeatFlux>().update(tst);
}

void TOPAZOcean::updateAfter(const TimestepTime& tst)
{
    mergeFluxes(tst);
    slabOcean.update(tst);
    sst = ModelArrayRef<Protected::SLAB_SST, RO>(getStore());
    sss = ModelArrayRef<Protected::SLAB_SSS, RO>(getStore());

    try {
        checkFields();
    } catch (const std::exception& e) {
        throw std::runtime_error("TOPAZOcean::updateAfter: " + std::string(e.what()));
    }
}

ModelState TOPAZOcean::getStatePrognostic() const
{
    ModelState state = IOceanBoundary::getStatePrognostic();
    return state.merge(slabOcean.getStatePrognostic());
}

ModelState TOPAZOcean::getStateDiagnostic() const
{
    ModelState state = IOceanBoundary::getStateDiagnostic();
    return state.merge(slabOcean.getStateDiagnostic());
}

void TOPAZOcean::setFilePath(const std::string& filePathIn) { filePath = filePathIn; }

void TOPAZOcean::setDirectory(const std::string& dir) { fileDirectory() = dir; }
const std::string& TOPAZOcean::getDirectory() { return fileDirectory(); }
const std::string TOPAZOcean::addDirectory(const std::string& file)
{
    return getDirectory() + "/" + file;
}

void TOPAZOcean::setData(const ModelState::DataMap& ms)
{
    IOceanBoundary::setData(ms);

    sstExt.resize();
    sssExt.resize();

    addChecks({
        { "sstExt", &sstExt },
        { "sssExt", &sssExt },
    });

    slabOcean.setData(ms);
}

void TOPAZOcean::updateTf(size_t i, const TimestepTime& tst) { tf[i] = (*pFreezingPoint)(sss[i]); }

const std::string topazFromNSName(const std::string& nsName)
{
    static const std::map<std::string, std::string> topazFromNS = {
            {"sst", "temperature"},
            {"sss", "salinity"},
            {"mld", "mlp"},
            {"u", "u"},
            {"v", "v"},
    };
    return topazFromNS.at(nsName);
}

std::string topazFilenameFromYearMonth(const std::string& topazName, size_t year, size_t month)
{
    static const std::vector<std::string> months = {
            "01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"};
    std::string depthString = (topazName == "u" || topazName == "v") ? "30m" : "3m";
    std::string filename = "TP4DAILY_" + std::to_string(year) + months[month-1] + "_" + depthString + ".nc";
    return TOPAZOcean::addDirectory(filename);
}

std::pair<double, double> topazFractionalIndices(double longitude, double latitude) {

    // x and y positions in the source grid for each point on the target grid
    const double dRad = radians(0.089828);
    const double lat0 = radians(90);
    const double lon0 = radians(-45.);
    const double topazNx = 760;
    const double topazNy = 1100;

    double cosLat = cos(radians(latitude));
    // Uses cos(lat0) = 0.
    double k = 2/dRad/(1+sin(radians(latitude)));
    double iFrac = k * cosLat * sin(radians(longitude) - lon0) + topazNx/2;
    double jFrac = -k * cosLat * cos(radians(longitude) - lon0) + topazNy/2;

    return {iFrac, jFrac};
}

ModelArray getTOPAZVarIndexData(const std::string& topazName, size_t year, size_t month, size_t day, const ModelArray& modelLon, const ModelArray& modelLat)
{
    std::string filePath = topazFilenameFromYearMonth(topazName, year, month);
    NetCDFForcings::Buffer tdata = NetCDFForcings::getFileIndexData(filePath, topazName, day-1);

    // x and y positions in the source grid for each point on the target grid
    double dRad = radians(0.09);
    double lat0 = radians(90);
    double lon0 = radians(0.);
    double topazNx = 761;
    double topazNy = 1101;

    ModelArray iFrac;
    ModelArray jFrac;
    ModelArray::MultiDim dim = modelLon.dimensions();
    for (int j = 0; j < dim[1]; ++j) {
        for (int i = 0; i < dim[0]; ++i) {
            auto [l, m] = topazFractionalIndices(modelLon(i, j), modelLat(i, j));
            iFrac(i, j) = l;
            jFrac(i, j) = m;
        }
    }
    return NetCDFForcings::maFromBuffer(tdata, iFrac, jFrac);
}

// Time interpolation happens here
ModelArray getTOPAZVarTimeData(const std::string& topazName, const TimePoint& time, const ModelArray& modelLon, const ModelArray& modelLat)
{
    ModelArray v1 = getTOPAZVarIndexData(topazName, time.year(), time.month(), time.day(), modelLon, modelLat);
    TimePoint t2 = time + Duration(3600);
    ModelArray v2 = getTOPAZVarIndexData(topazName, t2.year(), t2.month(), t2.day(), modelLon, modelLat);
    double f = t2.minute() / 60.;
    return v2 * f + v1 * (1-f);
}


const ModelArray TOPAZOcean::getData(const std::string& nsName, const TimePoint& time) const
{
    const ModelMetadata& meta = ModelMetadata::getInstance();
    return getTOPAZVarTimeData(topazFromNSName(nsName), time, meta.longitude(), meta.latitude());
}
} /* namespace Nextsim */
