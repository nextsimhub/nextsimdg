/*!
 * @file ModelComponent.hpp
 *
 * @date 1 Jul 2024
 * @author Tim Spain <timothy.spain@nersc.no>
 * @author Einar Ólason <einar.olason@nersc.no>
 */

#ifndef MODELCOMPONENT_HPP
#define MODELCOMPONENT_HPP

#include "include/Logged.hpp"
#include "include/MissingData.hpp"
#include "include/ModelArrayRef.hpp"
#include "include/ModelState.hpp"
#include "include/OutputSpec.hpp"
#include "include/TextTag.hpp"
#include "include/Time.hpp"

#include "ModelArrayRef.hpp"
#include "ModelArrayReferenceStore.hpp"
#include <functional>
#include <string>
#include <unordered_map>
#include <unordered_set>

namespace Nextsim {

class ModelComponent;

namespace Protected {
    // Prognostic model fields
    inline constexpr TextTag H_ICE = "H_ICE_cell"; // Ice thickness, cell average, m
    inline constexpr TextTag C_ICE = "C_ICE0"; // Ice concentration
    inline constexpr TextTag H_SNOW = "H_SNOW_cell"; // Snow depth, cell average, m
    inline constexpr TextTag T_ICE = "T_ICE0"; // Ice temperature, ˚C
    inline constexpr TextTag DAMAGE = "DAMAGE0"; // Ice damage 0–1
    // External data fields
    inline constexpr TextTag T_AIR = "T_AIR"; // Air temperature, ˚C
    inline constexpr TextTag DEW_2M = "DEW_2M"; // Dew point at 2 m, ˚C
    inline constexpr TextTag P_AIR = "P_AIR"; // sea level air pressure, Pa
    inline constexpr TextTag MIXRAT = "MIXRAT"; // water vapour mass mixing ratio
    inline constexpr TextTag SW_IN = "SW_IN"; // incoming shortwave flux, W m⁻²
    inline constexpr TextTag LW_IN = "LW_IN"; // incoming longwave flux, W m⁻²
    inline constexpr TextTag MLD = "MLD"; // mixed layer depth, m
    inline constexpr TextTag SNOW = "SNOWFALL"; // snow fall, kg m⁻² s⁻¹
    inline constexpr TextTag SSS = "SSS"; // sea surface salinity, PSU
    inline constexpr TextTag SST = "SST"; // sea surface temperature ˚C
    inline constexpr TextTag EXT_SSS
        = "EXT_SSS"; // sea surface salinity from coupling or forcing, PSU
    inline constexpr TextTag EXT_SST
        = "EXT_SST"; // sea surface temperature from coupling or forcing, ˚C
    inline constexpr TextTag EVAP_MINUS_PRECIP
        = "E-P"; // E-P atmospheric freshwater flux, kg s⁻¹ m⁻²
    // Derived fields, calculated once per timestep
    inline constexpr TextTag ML_BULK_CP = "CPML"; // Mixed layer bulk heat capacity J K⁻¹ m⁻²
    inline constexpr TextTag TF = "TF"; // Ocean freezing temperature, ˚C
    inline constexpr TextTag WIND_SPEED = "WIND_SPEED"; // Wind speed, m s⁻¹
    inline constexpr TextTag HTRUE_ICE = "HTRUE_ICE"; // Ice thickness, ice average, m
    inline constexpr TextTag HTRUE_SNOW = "HTRUE_SNOW"; // Snow thickness, ice average, m
    inline constexpr TextTag OCEAN_U = "OCEAN_U"; // x(east)-ward ocean current, m s⁻¹
    inline constexpr TextTag OCEAN_V = "OCEAN_V"; // y(north)-ward ocean current, m s⁻¹
    inline constexpr TextTag WIND_U = "WIND_U"; // x(east)-ward component of wind, m s⁻¹
    inline constexpr TextTag WIND_V = "WIND_V"; // y(north)-ward component of wind, m s⁻¹
    inline constexpr TextTag ICE_U = "ICE_U"; // x(east)-ward ice velocity, m s⁻¹
    inline constexpr TextTag ICE_V = "ICE_V"; // y(north)-ward ice velocity, m s⁻¹
    // Slab ocean fields
    inline constexpr TextTag SLAB_SST = "SLAB_SST"; // Slab ocean sea surface temperature, ˚C
    inline constexpr TextTag SLAB_SSS = "SLAB_SSS"; // Slab ocean sea surface salinity, ˚C
    inline constexpr TextTag SLAB_QDW
        = "SLAB_QDW"; // Slab ocean temperature nudging heat flux, W m⁻²
    inline constexpr TextTag SLAB_FDW
        = "SLAB_FDW"; // Slab ocean salinity nudging water flux, kg s⁻¹ m⁻²
}

namespace Shared {
    // Values of the prognostic fields updated during the timestep
    inline constexpr TextTag H_ICE = "H_ICE"; // Updated ice thickness, ice average, m
    inline constexpr TextTag C_ICE = "C_ICE"; // Updated ice concentration
    inline constexpr TextTag H_SNOW = "H_SNOW"; // Updated snow depth, ice average, m
    inline constexpr TextTag T_ICE = "T_ICE"; // Updated ice temperatures, ˚C
    inline constexpr TextTag DAMAGE = "DAMAGE"; // Updated damage 0–1
    // Heat fluxes
    inline constexpr TextTag Q_IA = "Q_IA"; // Ice to atmosphere heat flux W m⁻²
    inline constexpr TextTag Q_IC = "Q_IC"; // Ice conduction heat flux W m⁻²
    inline constexpr TextTag Q_IO = "Q_IO"; // Ice to ocean heat flux W m⁻²
    inline constexpr TextTag Q_OW = "Q_OW"; // Open water heat flux W m⁻²
    inline constexpr TextTag DQIA_DT
        = "DQIA_DT"; // Derivative of Qᵢₐ w.r.t. ice surface temperature  W m⁻² K⁻¹
    inline constexpr TextTag Q_PEN_SW = "Q_PEN_SW"; // Penetrating shortwave flux W m⁻²
    // Mass fluxes
    inline constexpr TextTag HSNOW_MELT = "HSNOW_MELT"; // Thickness of snow that melted, m
    // Atmospheric conditions
    inline constexpr TextTag SUBLIM = "SUBLIM"; // Upward sublimation rate kg m⁻² s⁻¹
    inline constexpr TextTag DELTA_HICE = "DELTA_HICE"; // Change in sea ice thickness, m
    inline constexpr TextTag DELTA_CICE = "DELTA_CICE"; // Change in sea ice concentration
    // Ice growth (that is not included above)
    inline constexpr TextTag NEW_ICE = "NEW_ICE"; // Volume of new ice formed [m]

}
/*!
 * A class encapsulating a component of the model. It also provide a method of
 * communicating data between ModelComponents using enums, static arrays of
 * pointers and the ModelArrayRef class.
 */
class ModelComponent {
public:
    typedef Logged::level OutputLevel;
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
     * @brief Returns the ModelArrayRef backing store.
     */
    static ModelArrayReferenceStore& getStore() { return store; }

protected:
    void registerModule();

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
    static std::set<std::string> advectedFields;

private:
    static ModelArrayReferenceStore store;
    static std::unordered_map<std::string, ModelComponent*> registeredModules;

    static size_t nOcean;
    static std::vector<size_t> oceanIndex;
};

} /* namespace Nextsim */

#endif /* MODELCOMPONENT_HPP */
