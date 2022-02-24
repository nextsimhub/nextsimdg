/*!
 * @file IIceOceanHeatFluxModule.hpp
 *
 * @date Feb 21, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef PHYSICS_SRC_MODULES_INCLUDE_IICEOCEANHEATFLUXMODULE_HPP
#define PHYSICS_SRC_MODULES_INCLUDE_IICEOCEANHEATFLUXMODULE_HPP

#include "include/Module.hpp"

#include "include/IIceOceanHeatFlux.hpp"

namespace Module {

template <> Module<Nextsim::IIceOceanHeatFlux>::map Module<Nextsim::IIceOceanHeatFlux>::functionMap;
class IIceOceanHeatFluxModule : public Module<Nextsim::IIceOceanHeatFlux> {
};

} /* namespace Module */

#endif /* PHYSICS_SRC_MODULES_INCLUDE_IICEOCEANHEATFLUXMODULE_HPP */
