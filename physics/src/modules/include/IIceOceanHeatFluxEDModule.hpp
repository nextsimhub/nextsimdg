/*!
 * @file IIceOceanHeatFluxModule.hpp
 *
 * @date Feb 21, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef PHYSICS_SRC_MODULES_INCLUDE_IICEOCEANHEATFLUXMODULE_HPP
#define PHYSICS_SRC_MODULES_INCLUDE_IICEOCEANHEATFLUXMODULE_HPP

#include "include/Module.hpp"

#include "include/IIceOceanHeatFluxED.hpp"

namespace Module {

template <> Module<Nextsim::IIceOceanHeatFluxED>::map Module<Nextsim::IIceOceanHeatFluxED>::functionMap;
class IIceOceanHeatFluxEDModule : public Module<Nextsim::IIceOceanHeatFluxED> {
};

} /* namespace Module */

#endif /* PHYSICS_SRC_MODULES_INCLUDE_IICEOCEANHEATFLUXMODULE_HPP */
