/*!
 * @file OceanStateModule.hpp
 *
 * @date May 11, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef OCEANSTATEMODULE_HPP
#define OCEANSTATEMODULE_HPP

#include "include/Module.hpp"

#include "include/OceanState.hpp"

namespace Module {

template <>
Module<Nextsim::OceanState>::map Module<Nextsim::OceanState>::functionMap;

template <>
Module<Nextsim::OceanState>::fn Module<Nextsim::OceanState>::spf;

class OceanStateModule : public Module<Nextsim::OceanState> {
    struct Constructor {
        Constructor();
    };
    static Constructor ctor;
};

} /* namespace Module */

#endif /* OCEANSTATEMODULE_HPP */
