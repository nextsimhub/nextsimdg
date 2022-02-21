/*!
 * @file IThermodynamicsModule.hpp
 *
 * @date Feb 21, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef PHYSICS_SRC_MODULES_INCLUDE_ITHERMODYNAMICSMODULE_HPP
#define PHYSICS_SRC_MODULES_INCLUDE_ITHERMODYNAMICSMODULE_HPP

#include "include/Module.hpp"

#include "include/IThermodynamics.hpp"

namespace Module {

class IThermodynamicsModule : public Module<Nextsim::IThermodynamics> {
};

} /* namespace Module */

#endif /* PHYSICS_SRC_MODULES_INCLUDE_ITHERMODYNAMICSMODULE_HPP */
