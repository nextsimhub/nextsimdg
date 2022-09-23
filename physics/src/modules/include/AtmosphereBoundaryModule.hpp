/*!
 * @file AtmosphereBoundaryModule.hpp
 *
 * @date Sep 23, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef ATMOSPHEREBOUNDARYMODULE_HPP
#define ATMOSPHEREBOUNDARYMODULE_HPP

#include "include/Module.hpp"

#include "include/IAtmosphereBoundary.hpp"

namespace Module {

template <> Module<Nextsim::IAtmosphereBoundary>::map Module<Nextsim::IAtmosphereBoundary>::functionMap;
class AtmosphereBoundaryModule : public Module<Nextsim::IAtmosphereBoundary> {
    struct Constructor {
    Constructor();
    };
    static Constructor ctor;
};

} /* namespace Module */

#endif /* ATMOSPHEREBOUNDARYMODULE_HPP */
