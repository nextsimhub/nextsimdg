/*!
 * @file ModelComponent.hpp
 *
 * @date Feb 28, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef MODELCOMPONENT_HPP
#define MODELCOMPONENT_HPP

#include "include/Logged.hpp"
#include "include/MissingData.hpp"
#include "include/ModelArrayRef.hpp"
#include "include/ModelState.hpp"
#include "include/OutputSpec.hpp"
#include "include/Time.hpp"

#include <functional>
#include <string>
#include <unordered_map>
#include <unordered_set>

namespace Nextsim {

class ModelComponent;

/*!
 * A class encapsulating a component of the model. It also provide a method of
 * communicating data between ModelComponents using enums, static arrays of
 * pointers and the ModelArrayRef class.
 */
class ModelComponent {
public:
    typedef Logged::level OutputLevel;

    enum class ProtectedArray {
        // Prognostic model fields
        H_ICE, // Ice thickness, cell average, m
        C_ICE, // Ice concentration
        H_SNOW, // Snow depth, cell average, m
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
        SSS, // sea surface salinity, PSU
        SST, // sea surface temperature, ˚C
        EXT_SSS, // sea surface salinity from coupling or forcing, PSU
        EXT_SST, // sea surface temperature from coupling or forcing, ˚C
        EVAP_MINUS_PRECIP, // E-P atmospheric freshwater flux, kg s⁻¹ m⁻²
        // Derived fields, calculated once per timestep
        ML_BULK_CP, // Mixed layer bulk heat capacity J K⁻¹ m⁻²
        TF, // Ocean freezing temperature, ˚C
        WIND_SPEED, // Wind speed, m s⁻¹
        WIND_U, // wind velocity x component, m s⁻¹
        WIND_V, // wind velocity y component, m s⁻¹
        HTRUE_ICE, // Ice thickness, ice average, m
        HTRUE_SNOW, // Snow thickness, ice average, m
        OCEAN_U, // x(east)-ward ocean current, m s⁻¹
        OCEAN_V, // y(north)-ward ocean current, m s⁻¹
        ICE_U, // x(east)-ward ice velocity, m s⁻¹
        ICE_V, // y(north)-ward ice velocity, m s⁻¹
        // Slab ocean fields
        SLAB_SST, // Slab ocean sea surface temperature, ˚C
        SLAB_SSS, // Slab ocean sea surface salinity, ˚C
        SLAB_QDW, // Slab ocean temperature nudging heat flux, W m⁻²
        SLAB_FDW, // Slab ocean salinity nudging water flux, kg s⁻¹ m⁻²
        COUNT // Count of enum values
    };

#ifdef DEBUG_MODELARRAYREF
    static const size_t SharedArrayOffset = static_cast<size_t>(ProtectedArray::COUNT);
#else
    static const size_t SharedArrayOffset = 0;
#endif
    enum class SharedArray {
        // Values of the prognostic fields updated during the timestep
        H_ICE = SharedArrayOffset, // Updated ice thickness, ice average, m
        C_ICE, // Updated ice concentration
        H_SNOW, // Updated snow depth, ice average, m
        T_ICE, // Updated ice temperatures, ˚C
        // Heat fluxes
        Q_IA, // Ice to atmosphere heat flux W m⁻²
        Q_IC, // Ice conduction heat flux W m⁻²
        Q_IO, // Ice to ocean heat flux W m⁻²
        Q_OW, // Open water heat flux W m⁻²
        DQIA_DT, // Derivative of Qᵢₐ w.r.t. ice surface temperature  W m⁻² K⁻¹
        Q_PEN_SW, // Short-wave flux penetrating the very surface of the ice W m⁻²
        // Mass fluxes
        HSNOW_MELT, // Thickness of snow that melted, m
        // Atmospheric conditions
        SUBLIM, // Upward sublimation rate kg m⁻² s⁻¹
        DELTA_HICE, // Change in sea ice thickness, m
        DELTA_CICE, // Change in sea ice concentration
        // Ice growth (that is not included above)
        NEW_ICE, // Volume of new ice formed [m]
        COUNT // Count of enum values
    };
    typedef std::function<void(size_t, const TimestepTime&)> IteratedFn;

    ModelComponent();
    virtual ~ModelComponent() = default;

    //! Returns the name of the component
    virtual std::string getName() const = 0;

    /*!
     * @brief Set the initial data of the component from the passed ModelState.
     *
     * @param state The ModelState containing the data to be set.
     */
    virtual void setData(const ModelState::DataMap& state) = 0;
    /*!
     * @brief Returns a ModelState from this component.
     *
     * @details The ModelState is map between field names and ModelArray data
     * arrays. The intention is to merge together different ModelSatates to
     * produce a combined state. The returned ModelState will include the
     * states of any subsidiary components held by the object. This is the
     * default level of output and should include all and only fields to be
     * placed in the restart file.
     */
    virtual ModelState getState() const = 0;
    /*!
     * @brief Returns a ModelState from this component at a specified level.
     *
     * @details See the zero argument version for more details. The output
     * levels reuse those defined in the Logged class. The default level is
     * NOTICE, and only levels such as INFO, DEBUG and TRACE should be used,
     * and should provide extra diagnostic fields.
     */
    virtual ModelState getState(const OutputLevel&) const = 0;

    /*!
     * @brief Returns the state of the ModelComponent and any ModelComponents
     * it depends on.
     *
     * @details Used to traverse the current tree of ModelComponents and return
     * the overall model state for diagnostic output.
     */
    virtual ModelState getStateRecursive(const OutputSpec& os) const
    {
        return os ? getState() : ModelState();
    }

    //! @brief Returns the names of all Type::H ModelArrays defined in this component.
    virtual std::unordered_set<std::string> hFields() const { return {}; }
    //! @brief Returns the names of all Type::U ModelArrays defined in this component.
    virtual std::unordered_set<std::string> uFields() const { return {}; }
    //! @brief Returns the names of all Type::V ModelArrays defined in this component.
    virtual std::unordered_set<std::string> vFields() const { return {}; }
    //! @brief Returns the names of all Type::Z ModelArrays defined in this component.
    virtual std::unordered_set<std::string> zFields() const { return {}; }

    static void setAllModuleData(const ModelState& stateIn);
    static ModelState getAllModuleState();
    static void unregisterAllModules();

    static void getAllFieldNames(std::unordered_set<std::string>& uF,
        std::unordered_set<std::string>& vF, std::unordered_set<std::string>& zF);

    /*!
     * @brief Registers a ModelArray into a SharedArray slot from outside any
     *        ModelComponent object. Intended for testing and debugging.
     */
    static void registerExternalSharedArray(SharedArray type, ModelArray* addr)
    {
        registerSharedArray(type, addr);
    }
    /*!
     * @brief Registers a ModelArray into a ProtectedArray slot from outside
     *        any ModelComponent object. Intended for testing and debugging.
     */
    static void registerExternalProtectedArray(ProtectedArray type, ModelArray* addr)
    {
        registerProtectedArray(type, addr);
    }

    /*!
     * @brief Returns a const reference to the store for SharedArray fields
     */
    static const MARBackingStore& getSharedArray() { return sharedArrays; }

    /*!
     * @brief Returns a const reference to the store for ProtectedArray fields
     */
    static const MARConstBackingStore& getProtectedArray() { return protectedArrays; }

protected:
    void registerModule();

    /*!
     * Adds a pointer to a slot into the SharedArray array.
     *
     * @param type The SharedArray slot to add the ModelArray into.
     * @param addr The address of the ModelArray to be shared.
     */
    static void registerSharedArray(SharedArray type, ModelArray* addr);

    /*!
     * Adds a pointer to a slot into the ProtectedArray array.
     *
     * @param type The ProtectedArray slot to add the ModelArray into.
     * @param addr The address of the ModelArray to be shared (read only).
     */
    static void registerProtectedArray(ProtectedArray type, const ModelArray* addr);

    inline static void overElements(IteratedFn fn, const TimestepTime& tst)
    {
        for (size_t i = 0; i < nOcean; ++i) {
            fn(oceanIndex[i], tst);
        }
    }

    /*!
     * @brief Sets the model-wide land-ocean mask (for HField arrays).
     * @param mask The HField ModelArray containing the mask data.
     *             0/false is land, >0 is sea.
     */
    static void setOceanMask(const ModelArray& mask);
    /*!
     * If there is no valid land mask, assume all points are ocean and
     * initialize accordingly.
     */
    static void noLandMask();

    /*!
     * @brief Returns a copy of the provided ModelArray, masked according to the
     * land-ocean mask.
     * @param data The data to be masked.
     */
    static ModelArray mask(const ModelArray& data);

    /*!
     * @brief Returns the ocean mask.
     */
    static const ModelArray& oceanMask();

protected:
    static ModelArray* p_oceanMaskH;

private:
    static MARBackingStore sharedArrays;
    static MARConstBackingStore protectedArrays;
    static std::unordered_map<std::string, ModelComponent*> registeredModules;

    static size_t nOcean;
    static std::vector<size_t> oceanIndex;
};

} /* namespace Nextsim */

#endif /* MODELCOMPONENT_HPP */
