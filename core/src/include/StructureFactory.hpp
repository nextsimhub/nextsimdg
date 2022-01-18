/*!
 * @file StructureFactory.hpp
 *
 * @date Jan 18, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef CORE_SRC_INCLUDE_STRUCTUREFACTORY_HPP
#define CORE_SRC_INCLUDE_STRUCTUREFACTORY_HPP

namespace Nextsim {

#include "include/ModuleLoader.hpp"
#include "include/IStructure.hpp"
#include "include/DevGrid.hpp"

class StructureFactory {
public:
    StructureFactory() = 0;
    ~StructureFactory() = 0;

    static std::shared_ptr<IStructure> generate(const std::string& structureName);
};

} /* namespace Nextsim */

#endif /* CORE_SRC_INCLUDE_STRUCTUREFACTORY_HPP */
