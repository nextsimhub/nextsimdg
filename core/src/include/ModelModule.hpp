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

    enum class ProtectedArray {
        // Prognostic model fields
        H_ICE, // Ice thickness, m
        C_ICE, // Ice concentration
        H_SNOW, // Snow depth, m
        T_ICE, // Ice temperature, ˚C
        // External data fields
        T_AIR, // Air temperature, ˚C
        DEW_2M, // Dew point at 2 m, ˚C
        P_AIR, // sea level air pressure, Pa
        MIXRAT, // water vapour mass mixing ratio
        SW_IN, // incoming shortwave flux, W m⁻²
        LW_IN, // incoming longwave flux, W m⁻²
        MLD, // mixed layer depth, m
        SNOW, // snow fall, kg m⁻² s⁻¹
    };
    enum class SharedArray {
        // Values of the prognostic fields updated during the timestep
        H_ICE, // Updated ice thickness, m
        C_ICE, // Updated ice concentration
        H_SNOW, // Updated snow depth, m
        T_ICE, // Updated ice temperatures, ˚C
        // Heat fluxes
        Q_IA, // Ice to atmosphere heat flux W m⁻²
        Q_IC, // Ice conduction heat flux W m⁻²
        Q_IO, // Ice to ocean heat flux W m⁻²
        Q_OW, // Open water heat flux W m⁻²
        DQIA_DT, // Derivative of Qᵢₐ w.r.t. ice surface temperature  W m⁻² K⁻¹
        // Atmospheric conditions
        SUBLIM, // Upward sublimation rate kg m⁻² s⁻¹
    };

    ModelModule();
    virtual ~ModelModule() = default;

    virtual std::string getName() const = 0;

    virtual void setData(const ModelState&) = 0;
    virtual ModelState getState() const = 0;
    virtual ModelState getState(const OutputLevel&) const = 0;

    virtual std::set<std::string> hFields() const { return {}; }
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
