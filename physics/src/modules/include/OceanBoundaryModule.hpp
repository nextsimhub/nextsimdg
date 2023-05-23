/*!
 * @file OceanBoundaryModule.hpp
 *
 * @date Sep 26, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef OCEANBOUNDARYMODULE_HPP
#define OCEANBOUNDARYMODULE_HPP

#include "include/Module.hpp"

#include "include/IOceanBoundary.hpp"

namespace Module {

template <> Module<Nextsim::IOceanBoundary>::map Module<Nextsim::IOceanBoundary>::functionMap;
class OceanBoundaryModule : public Module<Nextsim::IOceanBoundary> {
    struct Constructor {
        Constructor();
    };
    static Constructor ctor;
};

} /* namespace Module */

#endif /* OCEANBOUNDARYMODULE_HPP */
