/*!
 * @file DynamicsModule.hpp
 *
 * @date 6 Jan 2023
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef DYNAMICSMODULE_HPP
#define DYNAMICSMODULE_HPP

#include "include/Module.hpp"
#include "include/IDynamics.hpp"

namespace Module {

template <> Module<Nextsim::IDynamics>::map Module<Nextsim::IDynamics>::functionMap;
class DynamicsModule : public Module<Nextsim::IDynamics> {
    struct Constructor {
        Constructor();
    };
    static Constructor ctor;
};

} /* namespace Module */

#endif /* DYNAMICSMODULE_HPP */
