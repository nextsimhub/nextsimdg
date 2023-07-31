/*!
 * @file IIceOceanHeatFluxModule.hpp
 *
 * @date Feb 21, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef ICEOCEANHEATFLUXMODULE_HPP
#define ICEOCEANHEATFLUXMODULE_HPP

#include "include/Module.hpp"

#include "include/IIceOceanHeatFlux.hpp"

namespace Module {

template <> Module<Nextsim::IIceOceanHeatFlux>::map Module<Nextsim::IIceOceanHeatFlux>::functionMap;
class IceOceanHeatFluxModule : public Module<Nextsim::IIceOceanHeatFlux> {
    struct Constructor {
        Constructor();
    };
    static Constructor ctor;
};

} /* namespace Module */

#endif /* ICEOCEANHEATFLUXMODULE_HPP */
