/*!
 * @file VerticalIceGrowthModule.hpp
 *
 * @date Mar 17, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef VERTICALICEGROWTHMODULE_HPP
#define VERTICALICEGROWTHMODULE_HPP

#include "include/Module.hpp"

#include "include/IVerticalIceGrowth.hpp"

namespace Module {

template <>
Module<Nextsim::IVerticalIceGrowth>::map Module<Nextsim::IVerticalIceGrowth>::functionMap;
template <> Module<Nextsim::IVerticalIceGrowth>::fn Module<Nextsim::IVerticalIceGrowth>::spf;
class VerticalIceGrowthModule : public Module<Nextsim::IVerticalIceGrowth> {
    struct Constructor {
        Constructor();
    };
    static Constructor ctor;
};

} /* namespace Module */

#endif /* VERTICALICEGROWTHMODULE_HPP */
