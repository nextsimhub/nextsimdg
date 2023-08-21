/*!
 * @file IDiagnosticOutput.hpp
 *
 * @date May 25, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef IDIAGNOSTICOUTPUT_HPP
#define IDIAGNOSTICOUTPUT_HPP

#include "include/ModelComponent.hpp"
#include "include/ModelMetadata.hpp"
#include "include/ModelState.hpp"

#include <string>

namespace Nextsim {
class IDiagnosticOutput : public ModelComponent {
public:
    IDiagnosticOutput()
        : protectedArrayNames(generateProtectedNames())
        , sharedArrayNames(generateSharedNames())
        , protectedExternalNames(generateProtectedExternal())
        , sharedExternalNames(generateSharedExternal())
    {
    }
    virtual ~IDiagnosticOutput() = default;

    /*!
     * @brief Sets the output file name.
     *
     * @param fileName The file name to be set.
     */
    virtual void setFilenamePrefix(const std::string& filePrefix) = 0;

    /*!
     * @brief Outputs the passed ModelState.
     *
     * @param state The model state to be written out.
     * @param meta The model metadata for the the given state.
     */
    virtual void outputState(const ModelMetadata& meta) = 0;

    // Define some of the ModelComponent class functions
    // No data to be set
    void setData(const ModelState::DataMap& state) { }
protected:
    const std::map<std::string, ProtectedArray> protectedArrayNames;
    const std::map<std::string, SharedArray> sharedArrayNames;
    const std::map<std::string, std::string> protectedExternalNames;
    const std::map<std::string, std::string> sharedExternalNames;
private:
    /*
     * Using a pair of .ipp files to allow definitions of the externally visible
     * names of the fields to defined outside of an actual source file, even if the
     * definition file has a slightly odd format.
     */
    static inline std::map<std::string, std::string> generateProtectedExternal()
    {
        return {
#include "include/ProtectedArrayNames.ipp"
        };
    }
    static inline std::map<std::string, std::string> generateSharedExternal()
    {
        return {
#include "include/SharedArrayNames.ipp"
        };
    }

    // clang-format off
    // Using a macro protects against typos
    #define SHAREDEL(NAME) { #NAME, SharedArray::NAME }
    // Mapping from canonical SharedArray name to the SharedArray enum value
    static inline std::map<std::string, ModelComponent::SharedArray> generateSharedNames()
    {
        return {
            SHAREDEL(H_ICE), // Updated ice thickness, ice average, m
            SHAREDEL(C_ICE), // Updated ice concentration
            SHAREDEL(H_SNOW), // Updated snow depth, ice average, m
            SHAREDEL(T_ICE), // Updated ice temperatures, ˚C
            SHAREDEL(Q_IA), // Ice to atmosphere heat flux W m⁻²
            SHAREDEL(Q_IC), // Ice conduction heat flux W m⁻²
            SHAREDEL(Q_IO), // Ice to ocean heat flux W m⁻²
            SHAREDEL(Q_OW), // Open water heat flux W m⁻²
            SHAREDEL(DQIA_DT), // Derivative of Qᵢₐ w.r.t. ice surface temperature  W m⁻² K⁻¹
            SHAREDEL(HSNOW_MELT), // Thickness of snow that melted, m
            SHAREDEL(SUBLIM), // Upward sublimation rate kg m⁻² s⁻¹
            SHAREDEL(DELTA_HICE), // Change in sea ice thickness, m
            SHAREDEL(DELTA_CICE), // Change in sea ice concentration
            SHAREDEL(NEW_ICE), // Volume of new ice formed [m]
        };
    }
    #undef SHAREDEL

    // Using a macro protects against typos
    #define PROTEL(NAME) { #NAME, ProtectedArray::NAME }
    // Mapping from canonical ProtectedArray name to the ProtectedArray enum value
    static inline std::map<std::string, ModelComponent::ProtectedArray> generateProtectedNames()
    {
        return {
            PROTEL(H_ICE), // Ice thickness, cell average, m
            PROTEL(C_ICE), // Ice concentration
            PROTEL(H_SNOW), // Snow depth, cell average, m
            PROTEL(T_ICE), // Ice temperature, ˚C
            PROTEL(T_AIR), // Air temperature, ˚C
            PROTEL(DEW_2M), // Dew point at 2 m, ˚C
            PROTEL(P_AIR), // sea level air pressure, Pa
            PROTEL(MIXRAT), // water vapour mass mixing ratio
            PROTEL(SW_IN), // incoming shortwave flux, W m⁻²
            PROTEL(LW_IN), // incoming longwave flux, W m⁻²
            PROTEL(MLD), // mixed layer depth, m
            PROTEL(SNOW), // snow fall, kg m⁻² s⁻¹
            PROTEL(SSS), // sea surface salinity, PSU
            PROTEL(SST), // sea surface temperature ˚C
            PROTEL(EXT_SSS), // exterior sea surface salinity, PSU
            PROTEL(EXT_SST), // exterior sea surface temperature ˚C
            PROTEL(EVAP_MINUS_PRECIP), // E-P atmospheric freshwater flux, kg s⁻¹ m⁻²
            PROTEL(ML_BULK_CP), // Mixed layer bulk heat capacity J K⁻¹ m⁻²
            PROTEL(TF), // Ocean freezing temperature, ˚C
            PROTEL(WIND_SPEED), // Wind speed, m s⁻¹
            PROTEL(HTRUE_ICE), // Ice thickness, ice average, m
            PROTEL(HTRUE_SNOW), // Snow thickness, ice average, m
            PROTEL(OCEAN_U), // x(east)-ward ocean current, m s⁻¹
            PROTEL(OCEAN_V), // y(north)-ward ocean current, m s⁻¹
            PROTEL(SLAB_SSS), // slab ocean surface salinity, PSU
            PROTEL(SLAB_SST), // slab ocean surface temperature ˚C
            PROTEL(SLAB_QDW), // Slab ocean temperature nudging heat flux, W m⁻²
            PROTEL(SLAB_FDW), // Slab ocean salinity nudging water flux, kg s⁻¹ m⁻²
        };
    }
    #undef PROTEL
    // clang-format on

};
}
#endif /* IDIAGNOSTICOUTPUT_HPP */
