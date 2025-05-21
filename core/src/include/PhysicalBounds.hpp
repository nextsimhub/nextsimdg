/*!
 * @file PhysicalBounds.hpp
 *
 * @date 21 May 2025
 * @author Einar Ólason <einar.olason@nersc.no>
 */

#ifndef PHYSICALLIMITS_HPP
#define PHYSICALLIMITS_HPP

#include <map>

namespace Nextsim {

class PhysicalBounds;

/* List of reasonable physical bounds for all variables listed in ProtectedArrayNames.ipp and
 * SharedArrayNames.ipp */

// Shared arrays
class PhysicalBounds {
public:
    PhysicalBounds()
    {
        // Helper function to make the bounds.emplace code more readable
        auto addBounds = [this](const std::string& key, double min, double max) {
            bounds.emplace(key, std::make_pair(min, max));
        };

        // Protected arrays
        // Prognostic model fields
        addBounds(Protected::C_ICE, 0, 1); // Ice concentration
        addBounds(Protected::DAMAGE, 0, 1); // Ice thickness, cell average, m
        addBounds(Protected::H_ICE, 0, 50); // Ice thickness, cell average, m
        addBounds(Protected::H_SNOW, 0, 10); // Snow depth, cell average, m
        addBounds(Protected::ICE_U, -5, 5); // x(east)-ward ice velocity, m s⁻¹
        addBounds(Protected::ICE_V, -5, 5); // y(north)-ward ice velocity, m s⁻¹
        addBounds(Protected::SLAB_SSS, 0, 50); // Slab ocean surface salinity PSU
        addBounds(Protected::SLAB_SST, -5, 50); // Slab ocean surface temperature ˚C
        addBounds(Protected::SSS, 0, 50); // sea surface salinity, PSU
        addBounds(Protected::SST, -5, 50); // sea surface temperature ˚C
        addBounds(Protected::T_ICE, -100, 0); // Ice temperature, ˚C
        // External data fields
        addBounds(Protected::DEW_2M, -100, 80); // Dew point at 2 m, ˚C
        addBounds(
            Protected::EVAP_MINUS_PRECIP, -1, 1); // E-P atmospheric freshwater flux, kg s⁻¹ m⁻²
        addBounds(Protected::EXT_SSS, 0, 50); // External sea surface salinity PSU
        addBounds(Protected::EXT_SST, -5, 50); // External sea surface temperature ˚C
        addBounds(Protected::LW_IN, 0, 1e3); // incoming longwave flux, W m⁻²
        addBounds(Protected::MIXRAT, 0, 1); // water vapour mass mixing ratio
        addBounds(Protected::MLD, 1e-3, 12e3); // mixed layer depth, m
        addBounds(Protected::OCEAN_U, -1, 1); // x(east)-ward ocean current, m s⁻¹
        addBounds(Protected::OCEAN_V, -1, 1); // y(north)-ward ocean current, m s⁻¹
        addBounds(Protected::P_AIR, 700e2, 1200e2); // sea level air pressure, Pa
        addBounds(Protected::SNOW, 0, 1); // snow fall, kg m⁻² s⁻¹
        addBounds(Protected::SSH, -10, 10); // Slab ocean salinity nudging water flux, kg s⁻¹ m⁻²
        addBounds(Protected::SW_IN, -1e-6, 1e3); // incoming shortwave flux, W m⁻²
        addBounds(Protected::T_AIR, -100, 80); // Air temperature, ˚C
        addBounds(Protected::WIND_SPEED, 0, 200); // Wind speed, m s⁻¹
        addBounds(Protected::WIND_U, -200, 200); // wind velocity x component, m s⁻¹
        addBounds(Protected::WIND_V, -200, 200); // wind velocity y component, m s⁻¹
        // Derived fields, calculated once per timestep
        addBounds(Protected::HTRUE_ICE, 0, 50); // Ice thickness, ice average, m
        addBounds(Protected::HTRUE_SNOW, 0, 10); // Snow thickness, ice average, m
        addBounds(Protected::IO_STRESS_X, -10, 10); // Ice-ocean stress x(east) direction, Pa
        addBounds(Protected::IO_STRESS_Y, -10, 10); // Ice-ocean stress y(north) direction, Pa
        addBounds(Protected::ML_BULK_CP, 0, 1e11); // Mixed layer bulk heat capacity J K⁻¹ m⁻²
        addBounds(Protected::SLAB_FDW, -1, 1); // Slab ocean salinity nudging water flux, kg s⁻¹ m⁻²
        addBounds(
            Protected::SLAB_QDW, -1e5, 1e5); // Slab ocean temperature nudging heat flux, W m⁻²
        addBounds(Protected::TF, -5, 0); // Ocean freezing temperature, ˚C
        addBounds(Protected::T_SURF, -100, 0); // Ocean freezing temperature, ˚C

        // Shared arrays
        addBounds(Shared::C_ICE, 0, 1); // Updated ice concentration
        addBounds(Shared::C_ICE_DG, 0, 1); // Updated ice concentration
        addBounds(Shared::DAMAGE, 0, 1); // Updated ice thickness, ice average, m
        addBounds(Shared::DELTA_CICE, -0.1, 0.1); // Change in sea ice concentration
        addBounds(Shared::DELTA_HICE, -0.3, 0.3); // Change in sea ice thickness, m
        addBounds(Shared::DQIA_DT, -1e3,
            1e3); // Derivative of Qᵢₐ w.r.t. ice surface temperature  W m⁻² K⁻¹
        addBounds(Shared::HSNOW_MELT, 0, 0.1); // Thickness of snow that melted, m
        addBounds(Shared::H_ICE, 0, 50); // Updated ice thickness, ice average, m
        addBounds(Shared::H_ICE_DG, 0, 50); // Updated ice thickness, ice average, m
        addBounds(Shared::H_SNOW, 0, 10); // Updated snow depth, ice average, m
        addBounds(Shared::NEW_ICE, 0, 0.1); // Volume of new ice formed [m]
        addBounds(Shared::OW_STRESS_X, -10, 10); // x(east)-ward open ocean stress, Pa
        addBounds(Shared::OW_STRESS_Y, -10, 10); // y(east)-ward open ocean stress, Pa
        addBounds(Shared::Q_IA, -1e4, 1e4); // Ice to atmosphere heat flux W m⁻²
        addBounds(Shared::Q_IC, -1e4, 1e4); // Ice conduction heat flux W m⁻²
        addBounds(Shared::Q_IO, -1e8, 1e8); // Ice to ocean heat flux W m⁻²
        addBounds(Shared::Q_OW, -1e4, 1e4); // Open water heat flux W m⁻²
        addBounds(Shared::Q_PEN_SW, -1e-6, 100); // Penetrating shortwave flux W m⁻²
        addBounds(Shared::Q_SW_BASE, -10, 0); // Short-wave flux through ice base W m⁻²
        addBounds(Shared::Q_SW_OW, -1e3, 1e-6); // Short-wave flux into ice free ocean W m⁻²
        addBounds(Shared::SUBLIM, -1e-3, 1e-3); // Upward sublimation rate kg m⁻² s⁻¹
        addBounds(Shared::T_ICE, -100, 0); // Updated ice temperatures, ˚C
        addBounds(Shared::T_SURF, -100, 0); // Updated ice temperatures, ˚C
    }

    std::pair<double, double> getBounds(const std::string& key) const { return bounds.at(key); }

private:
    std::map<std::string, std::pair<double, double>> bounds;
};

}

#endif // PHYSICALBOUNDS_HPP
