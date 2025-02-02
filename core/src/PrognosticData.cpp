/*!
 * @file PrognosticData.cpp
 *
 * @date 02 Feb 2025
 * @author Tim Spain <timothy.spain@nersc.no>
 * @author Einar Ólason <einar.olason@nersc.no>
 */

#include "include/PrognosticData.hpp"

#include "include/Finalizer.hpp"
#include "include/ModelArrayRef.hpp"
#include "include/NextsimModule.hpp"
#include "include/gridNames.hpp"

namespace Nextsim {

const std::string PrognosticData::all = "ALL";

static const std::string pfx = "debug";
static const std::string fieldNamesKey = pfx + ".field_names";
static const std::string checkFieldsKey = pfx + ".check_fields";

static const std::map<int, std::string> keyMap = {
    { PrognosticData::FIELDNAMES_KEY, fieldNamesKey },
    { PrognosticData::CHECKFIELDS_KEY, checkFieldsKey },
};

static const std::string fieldNamesDefault = "hice,hsnow,cice,tice";
static const bool checkFieldsDefault = true;

PrognosticData::PrognosticData()
    : m_dt(1)
    , m_thick(ModelArray::Type::H)
    , m_conc(ModelArray::Type::H)
    , m_snow(ModelArray::Type::H)
    , m_tice(ModelArray::Type::Z)
    , m_damage(ModelArray::Type::H)
    , pAtmBdy(0)
    , pOcnBdy(0)
    , pDynamics(0)
    , bounds({
#include "include/PhysicalBounds.ipp"
      })
    , externalNames({
#include "include/ProtectedArrayNames.ipp"

#include "include/SharedArrayNames.ipp"

      })
{
    getStore().registerArray(Protected::H_ICE, &m_thick, RO);
    getStore().registerArray(Protected::C_ICE, &m_conc, RO);
    getStore().registerArray(Protected::H_SNOW, &m_snow, RO);
    getStore().registerArray(Protected::T_ICE, &m_tice, RO);
    getStore().registerArray(Protected::DAMAGE, &m_damage, RO);
}

void PrognosticData::configure()
{
    // Register finalizers before calling configure.
    Finalizer::registerUnique(Module::finalize<IAtmosphereBoundary>);
    Finalizer::registerUnique(Module::finalize<IOceanBoundary>);
    Finalizer::registerUnique(Module::finalize<IDynamics>);

    pAtmBdy = &Module::getImplementation<IAtmosphereBoundary>();
    tryConfigure(pAtmBdy);

    pOcnBdy = &Module::getImplementation<IOceanBoundary>();
    tryConfigure(pOcnBdy);

    pDynamics = &Module::getImplementation<IDynamics>();
    tryConfigure(pDynamics);

    tryConfigure(iceGrowth);
}

// Copies an HField from a source ModelArray that is either an HField or a DGField.
void copyMeanComponent(const ModelArray& source, ModelArray& sink)
{
    if (source.nComponents() > 1) {
        sink.setData(source.data().col(0));
    } else {
        sink = source;
    }
}

void PrognosticData::setData(const ModelState::DataMap& ms)
{

    if (ms.count(maskName)) {
        setOceanMask(ms.at(maskName));
    } else {
        noLandMask();
    }

    copyMeanComponent(ms.at(hiceName), m_thick);
    copyMeanComponent(ms.at(ciceName), m_conc);
    copyMeanComponent(ms.at(ticeName), m_tice);
    copyMeanComponent(ms.at(hsnowName), m_snow);
    // Damage is an optional field, and defaults to 1, if absent
    if (ms.count(damageName) > 0) {
        copyMeanComponent(ms.at(damageName), m_damage);
    } else {
        m_damage.resize();
        m_damage = 1.;
    }

    pAtmBdy->setData(ms);
    pOcnBdy->setData(ms);
    pDynamics->setData(ms);
    iceGrowth.setData(ms);

    if (getConfiguration(keyMap.at(CHECKFIELDS_KEY), checkFieldsDefault)) {
        setFieldsToCheck();
    }
}

// Go through the user supplied list of fields and get ModelArray* for each
void PrognosticData::setFieldsToCheck()
{
    // Populate a map of field name and pointers. Either through getStore.getAllData() or from a
    // user-supplied list
    std::unordered_map<std::string, const ModelArray*> storeData;
    if (const std::string listOfFields
        = getConfiguration(keyMap.at(FIELDNAMES_KEY), fieldNamesDefault);
        listOfFields == all) { // Check *all* the fields
        storeData = getStore().getAllData();
    } else { // Populate storeData with the fields listed
        std::istringstream fieldStream;
        fieldStream.str(listOfFields);
        for (std::string fieldName; std::getline(fieldStream, fieldName, ',');) {
            if (const ModelArray* fieldPointer
                = getStore().getArrayRef(externalNames.at(fieldName));
                fieldPointer == nullptr) {
                Logged::warning(
                    "PrognosticData: No field with the name \"" + fieldName + "\" was found.");
            } else {
                storeData.emplace(externalNames.at(fieldName), fieldPointer);
            }
        }
    }

    // Create a reverse look up for externalNames so the user can see the "user-friendly name"
    // in case of a crash
    std::map<std::string, std::string> reverseNames;
    for (const auto& ext : externalNames) {
        reverseNames[ext.second] = ext.first;
    }

    // Populate fieldsToCheck with user supplied fields
    for (const auto& x : storeData) {
        // Only check arrays that are in use (have been set/resized)
        if (x.second->data().rows() != 0) {
            fieldsToCheck.push_back({ reverseNames.at(x.first), x.second, bounds.at(x.first) });
        }
    }
}

void PrognosticData::update(const TimestepTime& tst)
{
    pOcnBdy->updateBefore(tst);
    pAtmBdy->update(tst);

    // Take the updated values of the true ice and snow thicknesses, and reset hice0 and hsnow0
    // IceGrowth updates its own fields during update
    iceGrowth.update(tst);
    updatePrognosticFields();

    pDynamics->update(tst);

    updatePrognosticFields();

    pOcnBdy->updateAfter(tst);

    checkFields(tst);
}

void PrognosticData::updatePrognosticFields()
{
    ModelArrayRef<Shared::H_ICE, RO> hiceTrueUpd(getStore());
    ModelArrayRef<Shared::C_ICE, RO> ciceUpd(getStore());
    ModelArrayRef<Shared::H_SNOW, RO> hsnowTrueUpd(getStore());
    ModelArrayRef<Shared::T_ICE, RO> ticeUpd(getStore());
    ModelArrayRef<Shared::DAMAGE, RO> damageUpd(getStore());

    // Calculate the cell average thicknesses
    HField hiceUpd = hiceTrueUpd * ciceUpd;
    HField hsnowUpd = hsnowTrueUpd * ciceUpd;

    m_thick.setData(hiceUpd);
    m_conc.setData(ciceUpd);
    m_snow.setData(hsnowUpd);
    m_tice.setData(ticeUpd);
    m_damage.setData(damageUpd);
}

// Gets all of the prognostic data, including that in the dynamics
ModelState PrognosticData::getState() const
{
    ModelArrayRef<Protected::SST> sst(getStore());
    ModelArrayRef<Protected::SSS> sss(getStore());

    // Get the prognostic data from the dynamics, including the full dynamics state
    ModelState dynamicsState = pDynamics->getState();
    // clang-format off
    ModelState localState = { {
                 { "mask", ModelArray(getOceanMask()) }, // make a copy
                 { "hice", mask(m_thick) },
                 { "cice", mask(m_conc) },
                 { "hsnow", mask(m_snow) },
                 { "tice", mask(m_tice) },
                 { "sst", mask(sst.data()) },
                 { "sss", mask(sss.data()) },
             },
        {} };
    // clang-format on

    // Use the dynamics values of any duplicated fields
    ModelState state(dynamicsState);
    state.merge(localState);

    return state;
}

// Recursively gets the data from all subcomponents
ModelState PrognosticData::getStateRecursive(const OutputSpec& os) const
{
    ModelState state;
    /* If allComponents is set on the OutputSpec, then for any duplicate fields, the subsystems
     * take priority, otherwise the fields held by PrognosticData itself. Note that std::map::merge
     * will not overwrite existing keys, so the first one that exists will survive.
     */
    if (os.allComponents()) {
        state.merge(pAtmBdy->getStateRecursive(os));
        state.merge(iceGrowth.getStateRecursive(os));
        state.merge(pDynamics->getStateRecursive(os));
        state.merge(getState());
    } else {
        state.merge(getState());
        state.merge(pAtmBdy->getStateRecursive(os));
        state.merge(iceGrowth.getStateRecursive(os));
        state.merge(pDynamics->getStateRecursive(os));
    }
    // OceanBoundary does not contribute to the output model state
    return os ? state : ModelState();
}

PrognosticData::HelpMap& PrognosticData::getHelpText(HelpMap& map, bool getAll)
{

    map["Debug"] = {
        { fieldNamesKey, ConfigType::STRING, {}, fieldNamesDefault, "",
            "Comma separated, space free list of fields to be checked if check_fields is true. "
            "The special value \""
                + all + "\" will output all available fields." },
        { checkFieldsKey, ConfigType::BOOLEAN, { "true", "false" },
            ConfigurationHelp::toString(checkFieldsDefault), "",
            "Switch to check if the fields listed in field_names are within a reasonable "
            "physical range and not NaN." },
    };
    return map;
}

PrognosticData::HelpMap& PrognosticData::getHelpRecursive(HelpMap& map, bool getAll)
{
    Module::getHelpRecursive<IAtmosphereBoundary>(map, getAll);
    Module::getHelpRecursive<IOceanBoundary>(map, getAll);
    Module::getHelpRecursive<IDynamics>(map, getAll);
    IceGrowth::getHelpRecursive(map, getAll);
    getHelpText(map, getAll);
    return map;
}

} /* namespace Nextsim */
