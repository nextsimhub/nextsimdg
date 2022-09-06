/*!
 * @file IIceAlbedoModule.hpp
 *
 * @date Feb 21, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef PHYSICS_SRC_MODULES_INCLUDE_IICEALBEDOMODULE_HPP
#define PHYSICS_SRC_MODULES_INCLUDE_IICEALBEDOMODULE_HPP

#include "include/ConfiguredModule.hpp"
#include "include/Module.hpp"

#include "include/IIceAlbedo.hpp"

namespace Module {

template <> Module<Nextsim::IIceAlbedo>::map Module<Nextsim::IIceAlbedo>::functionMap;
class IIceAlbedoModule : public Module<Nextsim::IIceAlbedo> {
    struct Constructor {
        Constructor();
    };
    static Constructor ctor;
};

} /* namespace Module */

#endif /* PHYSICS_SRC_MODULES_INCLUDE_IICEALBEDOMODULE_HPP */
