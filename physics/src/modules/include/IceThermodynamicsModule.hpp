/*!
 * @file IceThermodynamicsModule.hpp
 *
 * @date Mar 17, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef ICETHERMODYNAMICSMODULE_HPP
#define ICETHERMODYNAMICSMODULE_HPP

#include "include/Module.hpp"

#include "IIceThermodynamics.hpp"

namespace Module {

template <>
Module<Nextsim::IIceThermodynamics>::map Module<Nextsim::IIceThermodynamics>::functionMap;
template <> Module<Nextsim::IIceThermodynamics>::fn Module<Nextsim::IIceThermodynamics>::spf;
class IceThermodynamicsModule : public Module<Nextsim::IIceThermodynamics> {
    struct Constructor {
        Constructor();
    };
    static Constructor ctor;
};

} /* namespace Module */

#endif /* ICETHERMODYNAMICSMODULE_HPP */
