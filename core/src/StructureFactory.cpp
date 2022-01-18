/*!
 * @file StructureFactory.cpp
 *
 * @date Jan 18, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/StructureFactory.hpp"

#include <stdexcept>

namespace Nextsim {

std::shared_ptr<IStructure> StructureFactory::generate(const static std::string& structureName)
{
    ModuleLoader& loader = ModuleLoader::getLoader();
    std::string iStruct = "Nextsim::IStructure";
    std::shared_ptr<IStructure> shst;
    for (auto struc : loader.listImplementations(iStruct)) {
        loader.setImplementation(iStruct, struc);
        if (loader.getImplementation<IStructure>().structureType() == structureName) {
            shst = std::move(loader.getInstance<IStructure>());
        }
    }
    // If we reach here, then no valid handlers of the named structure were
    // found. Throw a invalid argument exception.
    std::string what = "StructureSelector: Invalid structure name (";
    what += structureName + ") provided.";
    throw std::invalid_argument(what);
}

} /* namespace Nextsim */
