/*!
 * @file IFreezingPointModule.hpp
 *
 * @date Feb 21, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef CORE_SRC_MODULES_INCLUDE_IFREEZINGPOINTMODULE_HPP
#define CORE_SRC_MODULES_INCLUDE_IFREEZINGPOINTMODULE_HPP

#include "include/Module.hpp"

#include "include/IFreezingPoint.hpp"

namespace Module {

template <> Module<Nextsim::IFreezingPoint>::map Module<Nextsim::IFreezingPoint>::functionMap;
class IFreezingPointModule : public Module<Nextsim::IFreezingPoint> {
    struct Constructor {
        Constructor();
    };
    static Constructor ctor;
};

} /* namespace Module */

#endif /* CORE_SRC_MODULES_INCLUDE_IFREEZINGPOINTMODULE_HPP */
