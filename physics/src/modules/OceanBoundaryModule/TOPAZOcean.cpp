/*!
 * @file TOPAZOcean.cpp
 *
 * @date 11 Feb 2025
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/TOPAZOcean.hpp"

#include "include/Finalizer.hpp"
#include "include/NextsimModule.hpp"
#include "include/ParaGridIO.hpp"
#include "include/constants.hpp"
#include "include/gridNames.hpp"

#include <include/PhysicalBounds.hpp>

namespace Nextsim {

std::string TOPAZOcean::filePath;
bool TOPAZOcean::doChecks;

static const std::string pfx = "TOPAZOcean";
static const std::string fileKey = pfx + ".file";
static const std::string checksKey = pfx + ".check_fields";

static const std::map<int, std::string> keyMap = {
    { TOPAZOcean::FILEPATH_KEY, fileKey },
    { TOPAZOcean::CHECKS_KEY, checksKey },
};

TOPAZOcean::TOPAZOcean()
    : sstExt(ModelArray::Type::H)
    , sssExt(ModelArray::Type::H)
{
}

ConfigurationHelp::HelpMap& TOPAZOcean::getHelpRecursive(HelpMap& map, bool getAll)
{
    map[pfx] = {
        { fileKey, ConfigType::STRING, {}, "", "",
            "Path to the processed NetCDF file providing the TOPAZ forcings." },
        { checksKey, ConfigType::BOOLEAN, {}, "false", "",
            "Check if the inputs are physically consistent." },
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
    doChecks = Configured::getConfiguration(keyMap.at(CHECKS_KEY), false);

    if (doChecks) {
        const PhysicalBounds bounds;
        // clang-format off
        fieldsToCheck.push_back({pfx+": sstExt", &sstExt, bounds.getBounds(Protected::EXT_SST)});
        fieldsToCheck.push_back({pfx+": sssExt", &sssExt, bounds.getBounds(Protected::EXT_SSS)});
        fieldsToCheck.push_back({pfx+": mld",    &mld,    bounds.getBounds(Protected::MLD)});
        fieldsToCheck.push_back({pfx+": u",      &u,      bounds.getBounds(Protected::OCEAN_U)});
        fieldsToCheck.push_back({pfx+": v",      &v,      bounds.getBounds(Protected::OCEAN_V)});
        fieldsToCheck.push_back({pfx+": ssh",    &ssh,    bounds.getBounds(Protected::SSH)});
        // clang-format on
    }

    slabOcean.configure();

    getStore().registerArray(Protected::EXT_SST, &sstExt, RO);
    getStore().registerArray(Protected::EXT_SSS, &sssExt, RO);
}

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

    cpml = Water::rho * Water::cp * mld;
    overElements(
        std::bind(&TOPAZOcean::updateTf, this, std::placeholders::_1, std::placeholders::_2),
        TimestepTime());

    pIOHeatFlux->update(tst);

    checkFields(tst);
}

void TOPAZOcean::updateAfter(const TimestepTime& tst)
{
    slabOcean.update(tst);
    sst = ModelArrayRef<Protected::SLAB_SST, RO>(getStore()).data();
    sss = ModelArrayRef<Protected::SLAB_SSS, RO>(getStore()).data();
}

void TOPAZOcean::setFilePath(const std::string& filePathIn) { filePath = filePathIn; }

void TOPAZOcean::setData(const ModelState::DataMap& ms)
{
    IOceanBoundary::setData(ms);

    sstExt.resize();
    sssExt.resize();
    slabOcean.setData(ms);
}

void TOPAZOcean::updateTf(size_t i, const TimestepTime& tst) { tf[i] = (*pFreezingPoint)(sss[i]); }

} /* namespace Nextsim */
