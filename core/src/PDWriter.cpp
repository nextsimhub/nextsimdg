/*!
 * @file PDWriter.cpp
 *
 * @date 14 Apr 2023
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/PrognosticData.hpp"

#include "include/ConfiguredModule.hpp"
#include "include/ModelMetadata.hpp"
#include "include/StructureFactory.hpp"

namespace Nextsim {
void PrognosticData::writeRestartFile(const std::string& filePath) const
{
    ConfigMap modelConfig;
    modelConfig.merge(getStateRecursive(true).config);
    modelConfig.merge(ConfiguredModule::getAllModuleConfigurations());
    ModelMetadata meta;
    meta.setConfig(modelConfig);
    StructureFactory::fileFromState(getState(), meta, filePath);
}
}
