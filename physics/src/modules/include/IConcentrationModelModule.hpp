/*!
 * @file IConcentrationModelModule.hpp
 *
 * @date Feb 21, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef PHYSICS_SRC_MODULES_INCLUDE_ICONCENTRATIONMODELMODULE_HPP
#define PHYSICS_SRC_MODULES_INCLUDE_ICONCENTRATIONMODELMODULE_HPP

#include "include/Module.hpp"

#include "include/IConcentrationModel.hpp"
namespace Module {

class IConcentrationModelModule : public Module<Nextsim::IConcentrationModel> {
};

} /* namespace Module */

#endif /* PHYSICS_SRC_MODULES_INCLUDE_ICONCENTRATIONMODELMODULE_HPP */
