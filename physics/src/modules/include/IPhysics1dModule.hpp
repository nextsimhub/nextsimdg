/*!
 * @file IPhysics1dModule.hpp
 *
 * @date Feb 18, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef PHYSICS_SRC_INCLUDE_IPHYSICS1DMODULE_HPP
#define PHYSICS_SRC_INCLUDE_IPHYSICS1DMODULE_HPP

#include "include/IPhysics1d.hpp"

#include "include/Module.hpp"

namespace Module {

class IPhysics1dModule : public Module<Nextsim::IPhysics1d> {
};

} /* namespace Module */

#endif /* PHYSICS_SRC_INCLUDE_IPHYSICS1DMODULE_HPP */
