/*!
 *
 * @author Tim Spain <timothy.spain@nersc.no>
*/

#include "include/PrognosticData.hpp"

#include "include/ConfiguredModule.hpp"
#include "include/ModelConfig.hpp"
#include "include/ModelMetadata.hpp"
#include "include/StructureFactory.hpp"

namespace Nextsim {
void PrognosticData::writeRestartFile(const std::string& filePath) const
{
    Logged::notice(std::string("  Writing state-based restart file: ") + filePath + '\n');

    ModelState state = getStatePrognostic();
    state.merge(ModelConfig::getConfig());

    auto& meta = ModelMetadata::getInstance();
    meta.affixCoordinates(state);

    StructureFactory::fileFromState(state, filePath, true);
}
}
