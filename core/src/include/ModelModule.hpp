/*!
 * @file ModelModule.hpp
 *
 * @date Feb 28, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef CORE_SRC_INCLUDE_MODELMODULE_HPP
#define CORE_SRC_INCLUDE_MODELMODULE_HPP

#include "include/Logged.hpp"
#include "include/ModelState.hpp"

#include <map>
#include <set>
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

    enum class ProtectedArray {
        H_ICE, // Ice thickness
        C_ICE, // Ice concentration
        H_SNOW, // Snow depth
        T_ICE, // Ice temperature
    };
    enum class SharedArray {
        Q_OW, // Open water heat flux
    };

    ModelModule();
    virtual ~ModelModule() = default;

    virtual std::string getName() const = 0;

    virtual void setData(const ModelState&) = 0;
    virtual ModelState getState() const = 0;
    virtual ModelState getState(const OutputLevel&) const = 0;

    virtual std::set<std::string> uFields() const { return {}; }
    virtual std::set<std::string> vFields() const { return {}; }
    virtual std::set<std::string> zFields() const { return {}; }

    static void setAllModuleData(const ModelState& stateIn);
    static ModelState getAllModuleState();
    static void unregisterAllModules();

    static void getAllFieldNames(
        std::set<std::string>& uF, std::set<std::string>& vF, std::set<std::string>& zF);

protected:
    void registerModule();

    static void registerSharedArray(SharedArray type, ModelArray* addr);
    static void requestSharedArray(SharedArray type, ModelArray** addr);

    static void registerProtectedArray(ProtectedArray type, const ModelArray* addr);
    static void requestProtectedArray(ProtectedArray, const ModelArray** addr);
private:
    static std::map<std::string, ModelModule*> registeredModules;
    static std::map<SharedArray, ModelArray*> registeredArrays;
    static std::map<SharedArray, std::set<ModelArray**>> reservedArrays;
    static std::map<ProtectedArray, const ModelArray*> registeredProtectedArrays;
    static std::map<ProtectedArray, std::set<const ModelArray**>> reservedProtectedArrays;
};

} /* namespace Nextsim */

#endif /* CORE_SRC_INCLUDE_MODELMODULE_HPP */
