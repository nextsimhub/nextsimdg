/*!
 * @file ModelMetadata.cpp
 *
 * @date Jun 29, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/ModelMetadata.hpp"

#include "include/StructureModule.hpp"

namespace Nextsim {

const std::string& ModelMetadata::structureName() const
{
    return Module::getImplementation<IStructure>().structureType();
}

} /* namespace Nextsim */
