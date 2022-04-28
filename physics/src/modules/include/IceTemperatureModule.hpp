/*!
 * @file IceTemperatureModule.hpp
 *
 * @date Apr 28, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef ICETEMPERATUREMODULE_HPP
#define ICETEMPERATUREMODULE_HPP

#include "include/Module.hpp"

#include "include/IIceTemperature.hpp"

namespace Module {

template <> Module<Nextsim::IIceTemperature>::map Module<Nextsim::IIceTemperature>::functionMap;
template <> Module<Nextsim::IIceTemperature>::fn Module<Nextsim::IIceTemperature>::spf;
class IceTemperatureModule : public Module<Nextsim::IIceTemperature> {
    struct Constructor {
        Constructor();
    };
    static Constructor ctor;
};

} /* namespace Module */

#endif /* ICETEMPERATUREMODULE_HPP */
