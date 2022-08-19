/*!
 * @file OceanState.cpp
 *
 * @date May 9, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/OceanState.hpp"

#include "include/constants.hpp"

namespace Nextsim {

static const std::map<ModelComponent::ProtectedArray, std::string> fieldNames = {
    { ModelComponent::ProtectedArray::SST, "Sea surface temperature" },
    { ModelComponent::ProtectedArray::SSS, "Sea surface salinity" },
    { ModelComponent::ProtectedArray::MLD, "Mixed layer depth" },
    { ModelComponent::ProtectedArray::TF, "Sea surface freezing temperature" },
    { ModelComponent::ProtectedArray::ML_BULK_CP, "Mixed layer bulk heat capacity" },
};

OceanState::OceanState()
    : tfImpl(nullptr)
    , sst(ModelArray::Type::H)
    , sss(ModelArray::Type::H)
    , mld(ModelArray::Type::H)
    , tf(ModelArray::Type::H)
    , cpml(ModelArray::Type::H)
{
    registerProtectedArray(ProtectedArray::SST, &sst);
    registerProtectedArray(ProtectedArray::SSS, &sss);
    registerProtectedArray(ProtectedArray::MLD, &mld);

    registerProtectedArray(ProtectedArray::TF, &tf);
    registerProtectedArray(ProtectedArray::ML_BULK_CP, &cpml);
}

void OceanState::setData(const ModelState::DataMap&) { }

ModelState OceanState::getState() const
{
    return {{
        { fieldNames.at(ProtectedArray::SST), sst },
        { fieldNames.at(ProtectedArray::SSS), sss },
        { fieldNames.at(ProtectedArray::MLD), mld },
        { fieldNames.at(ProtectedArray::TF), tf },
        { fieldNames.at(ProtectedArray::ML_BULK_CP), cpml },
    }, {}};
}

ModelState OceanState::getState(const OutputLevel& lvl) const { return getState(); }

std::string OceanState::getName() const { return "OceanState"; }

std::unordered_set<std::string> OceanState::hFields() const
{
    std::unordered_set<std::string> fields;
    for (const auto typeName : fieldNames) {
        fields.insert(typeName.second);
    }
    return fields;
}

void OceanState::configure()
{
    tfImpl = &Module::getImplementation<IFreezingPoint>();
    tryConfigure(tfImpl);
}

OceanState::HelpMap& OceanState::getHelpRecursive(HelpMap& map, bool getAll)
{
    return Module::getHelpRecursive<IFreezingPoint>(map, getAll);
}

void OceanState::update(const TimestepTime& tst)
{
    // Mixed layer heat capacity per unit area
    cpml = mld * Water::rho * Water::cp;
    // Sea surface freezing point
    tf.resize();
    overElements(std::bind(&OceanState::updateFreezingPoint, this, std::placeholders::_1,
                     std::placeholders::_2),
        tst);
    // Derived class updates
    updateSpecial(tst);
}

void OceanState::updateFreezingPoint(size_t i, const TimestepTime&) { tf[i] = (*tfImpl)(sss[i]); }

} /* namespace Nextsim */
