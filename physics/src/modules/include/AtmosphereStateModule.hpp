/*!
 * @file AtmosphereStateModule.hpp
 *
 * @date May 11, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef ATMOSPHERESTATEMODULE_HPP
#define ATMOSPHERESTATEMODULE_HPP

#include "include/Module.hpp"

#include "include/AtmosphereState.hpp"

namespace Module {

template <>
Module<Nextsim::AtmosphereState>::map Module<Nextsim::AtmosphereState>::functionMap;

template <>
Module<Nextsim::AtmosphereState>::fn Module<Nextsim::AtmosphereState>::spf;

class AtmosphereStateModule : public Module<Nextsim::AtmosphereState> {
    struct Constructor {
        Constructor();
    };
    static Constructor ctor;
};

} /* namespace Module */

#endif /* ATMOSPHERESTATEMODULE_HPP */
