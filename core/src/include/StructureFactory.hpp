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
    static std::shared_ptr<IStructure> generate(const std::string& structureName);

private:
    StructureFactory() = default;
};

} /* namespace Nextsim */

#endif /* CORE_SRC_INCLUDE_STRUCTUREFACTORY_HPP */
