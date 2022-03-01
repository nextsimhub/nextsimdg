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

#include <map>
#include <string>

namespace Nextsim {

class ModelModule;

class ModelModule {
public:
    typedef Logged::level OutputLevel;
    typedef ModelArray HField;
    typedef ModelArray UField;
    typedef ModelArray VField;
    typedef ModelArray ZField; // This needs to be made into a 3d field

    ModelModule();
    virtual ~ModelModule() = default;

    virtual std::string getName() const = 0;

    virtual void setData(const ModelState&) = 0;
    virtual ModelState getState() const = 0;
    virtual ModelState getState(const OutputLevel&) const = 0;

    static void setAllModuleData(const ModelState& stateIn);
    static ModelState getAllModuleState();
    static void unregisterAllModules();

protected:
    void registerModule();

private:
    static std::map<std::string, ModelModule*> registeredModules;
};

} /* namespace Nextsim */

#endif /* CORE_SRC_INCLUDE_MODELMODULE_HPP_ */
