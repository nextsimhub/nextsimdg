/*!
 * @file ModelModule.hpp
 *
 * @date Feb 28, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef CORE_SRC_INCLUDE_MODELMODULE_HPP_
#define CORE_SRC_INCLUDE_MODELMODULE_HPP_

#include "include/Logged.hpp"

namespace Nextsim {

class ModelState;

class ModelModule {
public:
    ModelModule();
    virtual ~ModelModule() = default;

    void setData(ModelState&);
    ModelState getState();
    ModelState getState(Logged::level&);
};

} /* namespace Nextsim */

#endif /* CORE_SRC_INCLUDE_MODELMODULE_HPP_ */
