/*!
 * @file PDWriter.cpp
 *
 * @date 14 Apr 2023
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/PrognosticData.hpp"

#include "include/ConfiguredModule.hpp"
#include "include/ModelConfig.hpp"
#include "include/ModelMetadata.hpp"
#include "include/StructureFactory.hpp"

namespace Nextsim {
void PrognosticData::writeRestartFile(
    const std::string& filePath, const ModelMetadata& metadata) const
{
    Logged::notice(std::string("  Writing state-based restart file: ") + filePath + '\n');

    ModelState state = getStatePrognostic();
    state.merge(ModelConfig::getConfig());

    ModelMetadata meta(metadata);
    meta.affixCoordinates(state);

    StructureFactory::fileFromState(state, meta, filePath, true);
}
}
