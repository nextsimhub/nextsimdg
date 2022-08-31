/*!
 * @file AtmosphereState.cpp
 *
 * @date May 9, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/AtmosphereState.hpp"

#include <vector>
namespace Nextsim {

static const std::map<ModelComponent::ProtectedArray, std::string> fieldNames = {
    { ModelComponent::ProtectedArray::T_AIR, "Air temperature" },
    { ModelComponent::ProtectedArray::DEW_2M, "Dew point temperature at 2 m" },
    { ModelComponent::ProtectedArray::P_AIR, "Air pressure" },
    { ModelComponent::ProtectedArray::MIXRAT, "Water vapour mixing ratio" },
    { ModelComponent::ProtectedArray::SW_IN, "Incident shortwave flux" },
    { ModelComponent::ProtectedArray::LW_IN, "Incident longwave flux" },
    { ModelComponent::ProtectedArray::SNOW, "Snowfall rate" },
    { ModelComponent::ProtectedArray::WIND_SPEED, "Wind speed at 10 m" },
};

AtmosphereState::AtmosphereState()
{
    registerProtectedArray(ProtectedArray::T_AIR, &tair);
    registerProtectedArray(ProtectedArray::DEW_2M, &tdew);
    registerProtectedArray(ProtectedArray::P_AIR, &pair);
    registerProtectedArray(ProtectedArray::MIXRAT, &rmix);
    registerProtectedArray(ProtectedArray::SW_IN, &sw_in);
    registerProtectedArray(ProtectedArray::LW_IN, &lw_in);
    registerProtectedArray(ProtectedArray::SNOW, &snowfall);

    registerProtectedArray(ProtectedArray::WIND_SPEED, &windSpeed);
}

void AtmosphereState::setData(const ModelState::DataMap&) { }

ModelState AtmosphereState::getState() const
{
    return { {
                 { fieldNames.at(ProtectedArray::T_AIR), tair },
                 { fieldNames.at(ProtectedArray::DEW_2M), tdew },
                 { fieldNames.at(ProtectedArray::P_AIR), pair },
                 { fieldNames.at(ProtectedArray::MIXRAT), rmix },
                 { fieldNames.at(ProtectedArray::SW_IN), sw_in },
                 { fieldNames.at(ProtectedArray::LW_IN), lw_in },
                 { fieldNames.at(ProtectedArray::SNOW), snowfall },
                 { fieldNames.at(ProtectedArray::WIND_SPEED), windSpeed },
             },
        {} };
}

ModelState AtmosphereState::getState(const OutputLevel& lvl) const { return getState(); }
ModelState AtmosphereState::getStateRecursive(const OutputSpec& os) const
{
    ModelState state(getState());
    return os ? state : ModelState();
}

std::string AtmosphereState::getName() const { return "AtmosphereState"; }

std::unordered_set<std::string> AtmosphereState::hFields() const
{
    std::unordered_set<std::string> fields;
    for (const auto typeName : fieldNames) {
        fields.insert(typeName.second);
    }
    return fields;
}

void AtmosphereState::update(const TimestepTime& tst)
{
    // Calculate wind speed from velocity components
    updateSpecial(tst);
}

} /* namespace Nextsim */
