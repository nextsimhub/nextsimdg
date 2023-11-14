/*!
 * @file IceAlbedoModule.hpp
 *
 * @date Feb 21, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef ICEALBEDOMODULE_HPP
#define ICEALBEDOMODULE_HPP

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

#endif /* ICEALBEDOMODULE_HPP */
