/*!
 * @file LateralIceSpreadModule.hpp
 *
 * @date Apr 8, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef LATERALICESPREADMODULE_HPP
#define LATERALICESPREADMODULE_HPP

#include "include/Module.hpp"

#include "include/ILateralIceSpread.hpp"

namespace Module {

template <> Module<Nextsim::ILateralIceSpread>::map Module<Nextsim::ILateralIceSpread>::functionMap;
class LateralIceSpreadModule : public Module<Nextsim::ILateralIceSpread> {
    struct Constructor {
        Constructor();
    };
    static Constructor ctor;
};

} /* namespace Module */

#endif /* LATERALICESPREADMODULE_HPP */
