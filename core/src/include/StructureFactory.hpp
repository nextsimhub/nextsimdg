/*!
 * @file StructureFactory.hpp
 *
 * @date Jan 18, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef CORE_SRC_INCLUDE_STRUCTUREFACTORY_HPP
#define CORE_SRC_INCLUDE_STRUCTUREFACTORY_HPP

#include "include/ModuleLoader.hpp"
#include "include/IStructure.hpp"
#include "include/DevGrid.hpp"

namespace Nextsim {

class StructureFactory {
public:
    /*!
     * @brief Returns a shared_ptr to a instance of IStructure which matches
     * the passed structure name.
     *
     * @param structureName the name of the structure that should provide the
     *        implementation.
     */
    static std::shared_ptr<IStructure> generate(const std::string& structureName);

    /*!
     * @brief Returns a shared_ptr to a instance of IStructure which implements
     * the structure named in the passed NetCDF file.
     *
     * @param filePath the name of the file to be read.
     */
    static std::shared_ptr<IStructure> generateFromFile(const std::string& fielPath);
private:
    StructureFactory() = default;
};

} /* namespace Nextsim */

#endif /* CORE_SRC_INCLUDE_STRUCTUREFACTORY_HPP */
