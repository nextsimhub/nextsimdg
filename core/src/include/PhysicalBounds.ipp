/*!
 * @file PrognosticData.cpp
 *
 * @date 30 Jan 2025
 * @author Einar Ólason <einar.olason@nersc.no>
 */

/* List of reasonable physical bounds for all variables listed in ProtectedArrayNames.ipp and
 * SharedArrayNames.ipp */

// Shared arrays
{ "DAMAGE0", {0, 1} }, // Ice thickness, cell average, m
{ "H_ICE_cell", {0, 50} }, // Ice thickness, cell average, m
{ "C_ICE0", {0, 1} }, // Ice concentration
{ "H_SNOW_cell", {0, 10} }, // Snow depth, cell average, m
{ "T_ICE0", {-100, 0} }, // Ice temperature, ˚C
{ "T_AIR", {-100, 80} }, // Air temperature, ˚C
{ "DEW_2M", {-100, 80} }, // Dew point at 2 m, ˚C
{ "P_AIR", {700e2, 1200e2} }, // sea level air pressure, Pa
{ "MIXRAT", {0, 1} }, // water vapour mass mixing ratio
{ "SW_IN", {-1e-6, 1e3} }, // incoming shortwave flux, W m⁻²
{ "LW_IN", {0, 1e3} }, // incoming longwave flux, W m⁻²
{ "MLD", {1e-3, 12e3} }, // mixed layer depth, m
{ "SNOWFALL", {0, 1} }, // snow fall, kg m⁻² s⁻¹
{ "SSS", {0, 50} }, // sea surface salinity, PSU
{ "SST", {-5, 50} }, // sea surface temperature ˚C
{ "EXT_SST", {-5, 50} }, // External sea surface temperature ˚C
{ "EXT_SSS", {0, 50} }, // External sea surface salinity PSU
{ "E-P", {-1, 1} }, // E-P atmospheric freshwater flux, kg s⁻¹ m⁻²
{ "CPML", {0,  1e11} }, // Mixed layer bulk heat capacity J K⁻¹ m⁻²
{ "TF", {-5, 0} }, // Ocean freezing temperature, ˚C
{ "WIND_SPEED", {0, 150} }, // Wind speed, m s⁻¹
{ "WIND_U", {-150, 150} }, // wind velocity x component, m s⁻¹
{ "WIND_V", {-150, 150} }, // wind velocity y component, m s⁻¹
{ "HTRUE_ICE", {0, 50} }, // Ice thickness, ice average, m
{ "HTRUE_SNOW", {0, 10} }, // Snow thickness, ice average, m
{ "OCEAN_U", {-1, 1} }, // x(east)-ward ocean current, m s⁻¹
{ "OCEAN_V", {-1, 1} }, // y(north)-ward ocean current, m s⁻¹
{ "ICE_U", {-5, 5} }, // x(east)-ward ice velocity, m s⁻¹
{ "ICE_V", {-5, 5} }, // y(north)-ward ice velocity, m s⁻¹
{ "SLAB_SST", {-5, 50} }, // Slab ocean surface temperature ˚C
{ "SLAB_SSS", {0, 50} }, // Slab ocean surface salinity PSU
{ "SLAB_QDW", {-1e5, 1e5} }, // Slab ocean temperature nudging heat flux, W m⁻²
{ "SLAB_FDW", {-1, 1} }, // Slab ocean salinity nudging water flux, kg s⁻¹ m⁻²
{ "SSH", {-10, 10} }, // Slab ocean salinity nudging water flux, kg s⁻¹ m⁻²
{ "IO_STRESS_X", {-10, 10} }, // Ice-ocean stress x(east) direction, Pa
{ "IO_STRESS_Y", {-10, 10} }, // Ice-ocean stress y(north) direction, Pa

// Protected arrays
{ "DAMAGE", {0, 1} }, // Updated ice thickness, ice average, m
{ "H_ICE", {0, 50} }, // Updated ice thickness, ice average, m
{ "C_ICE", {0, 1} }, // Updated ice concentration
{ "H_SNOW", {0, 10} }, // Updated snow depth, ice average, m
{ "T_ICE", {-100, 0} }, // Updated ice temperatures, ˚C
{ "Q_IA", {-1e4, 1e4} }, // Ice to atmosphere heat flux W m⁻²
{ "Q_IC", {-1e4, 1e4} }, // Ice conduction heat flux W m⁻²
{ "Q_IO", {-1e8, 1e8} }, // Ice to ocean heat flux W m⁻²
{ "Q_OW", {-1e4, 1e4} }, // Open water heat flux W m⁻²
{ "DQIA_DT", {-1e3, 1e3} }, // Derivative of Qᵢₐ w.r.t. ice surface temperature  W m⁻² K⁻¹
{ "Q_PEN_SW", {-1e-6, 100} }, // Penetrating shortwave flux W m⁻²
{ "HSNOW_MELT", {0, 0.1} }, // Thickness of snow that melted, m
{ "SUBLIM", {-1e-3, 1e-3} }, // Upward sublimation rate kg m⁻² s⁻¹
{ "DELTA_HICE", {-0.1, 0.1} }, // Change in sea ice thickness, m
{ "DELTA_CICE", {-0.1, 0.1} }, // Change in sea ice concentration
{ "NEW_ICE", {0, 0.1} }, // Volume of new ice formed [m]


