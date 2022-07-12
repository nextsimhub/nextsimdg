/*!
 * @file FluxCalculationModule.hpp
 *
 * @date Apr 29, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef FLUXCALCULATIONMODULE_HPP
#define FLUXCALCULATIONMODULE_HPP

#include "include/Module.hpp"

#include "include/IFluxCalculation.hpp"
namespace Module {

template <> Module<Nextsim::IFluxCalculation>::map Module<Nextsim::IFluxCalculation>::functionMap;

template <> Module<Nextsim::IFluxCalculation>::fn Module<Nextsim::IFluxCalculation>::spf;

class FluxCalculationModule : public Module<Nextsim::IFluxCalculation> {
    struct Constructor {
        Constructor();
    };
    static Constructor ctor;
};

} /* namespace Module */

#endif /* FLUXCALCULATIONMODULE_HPP */
