/*!
 * @file ModelModule.hpp
 *
 * @date Feb 28, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef CORE_SRC_INCLUDE_MODELMODULE_HPP_
#define CORE_SRC_INCLUDE_MODELMODULE_HPP_

#include "include/Logged.hpp"
#include "include/ModelState.hpp"

#include <functional>
#include <map>
#include <string>

namespace Nextsim {

class ModelModule;

class ModelModule {
public:
    ModelModule();
    virtual ~ModelModule() = default;

    virtual std::string getName() const = 0;

    virtual void setData(const ModelState&) = 0;
    virtual ModelState getState() const = 0;
    virtual ModelState getState(Logged::level&) const = 0;

    static void setAllModuleData(const ModelState& stateIn);
    static ModelState getAllModuleState();

protected:
    void registerModule();

private:
    static std::map<std::string, std::reference_wrapper<ModelModule>> registeredModules;
};

} /* namespace Nextsim */

#endif /* CORE_SRC_INCLUDE_MODELMODULE_HPP_ */
