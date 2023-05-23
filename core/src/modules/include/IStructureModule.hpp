/*!
 * @file IStructureModule.hpp
 *
 * @date Feb 18, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef CORE_SRC_MODULES_INCLUDE_ISTRUCTUREMODULE_HPP
#define CORE_SRC_MODULES_INCLUDE_ISTRUCTUREMODULE_HPP

#include "include/Module.hpp"

#include "include/IStructure.hpp"

namespace Module {

template <> Module<Nextsim::IStructure>::map Module<Nextsim::IStructure>::functionMap;
class IStructureModule : public Module<Nextsim::IStructure> {
};

} /* namespace Module */

#endif /* CORE_SRC_MODULES_INCLUDE_ISTRUCTUREMODULE_HPP */
