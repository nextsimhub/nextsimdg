/*!
 * @file StructureFactory.hpp
 *
 * @date Jan 18, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef STRUCTUREFACTORY_HPP
#define STRUCTUREFACTORY_HPP

#include "include/IStructure.hpp"

#include "include/ModelMetadata.hpp"
#include "include/ModelState.hpp"

#include <string>

namespace Nextsim {

class StructureFactory {
public:
    /*!
     * @brief Returns the ModelState of the named restart file.
     *
     * @param filePath the name of the file to be read.
     */
    static ModelState stateFromFile(const std::string& filePath);

    /*!
     * @brief Takes a ModelState and a template file name to write the state
     * out to a target file path.
     *
     * @param state the ModelState to be written.
     * @param filePath the path for the file to be written to.
     */
    static void fileFromState(const ModelState& state, const ModelMetadata& meta,
        const std::string& filePath, bool isRestart = false);

    static void finaliseAllFiles();

private:
    StructureFactory() = default;
};

} /* namespace Nextsim */

#endif /* CORE_SRC_INCLUDE_STRUCTUREFACTORY_HPP */
