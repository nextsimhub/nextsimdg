/*!
 * @file ProtectedArrayNames.ipp
 *
 * @date 23 Aug 2024
 * @author Tim Spain <timothy.spain@nersc.no>
 * @author Einar Ólason <einar.olason@nersc.no>
 */

// External ProtectedArray names must be unique wrt to the external SharedArray names as well
{ "dew2m", "DEW_2M" }, // Dew point at 2 m, ˚C
    { "eminusp", "E-P" }, // E-P atmospheric freshwater flux, kg s⁻¹ m⁻²
    { "fdw", "SLAB_FDW" }, // Slab ocean salinity nudging water flux, kg s⁻¹ m⁻²
    { "hice_true_pro", "HTRUE_ICE" }, // Ice thickness, ice average, m
    { "hsnow_true_pro", "HTRUE_SNOW" }, // Snow thickness, ice average, m
    { "lw_in", "LW_IN" }, // incoming longwave flux, W m⁻²
    { "mixrat", "MIXRAT" }, // water vapour mass mixing ratio
    { "mlcp", "CPML" }, // Mixed layer bulk heat capacity J K⁻¹ m⁻²
    { "mld", "MLD" }, // mixed layer depth, m
    { "ocean_u", "OCEAN_U" }, // x(east)-ward ocean current, m s⁻¹
    { "ocean_v", "OCEAN_V" }, // y(north)-ward ocean current, m s⁻¹
    { "pair", "P_AIR" }, // sea level air pressure, Pa
    { "qdw", "SLAB_QDW" }, // Slab ocean temperature nudging heat flux, W m⁻²
    { "snowfall", "SNOWFALL" }, // snow fall, kg m⁻² s⁻¹
    { "ssh", "SSH" }, // Slab ocean salinity nudging water flux, kg s⁻¹ m⁻²
    { "sss", "SSS" }, // sea surface salinity, PSU
    { "sss_ext", "EXT_SSS" }, // External sea surface salinity PSU
    { "sst", "SST" }, // sea surface temperature ˚C
    { "sst_ext", "EXT_SST" }, // External sea surface temperature ˚C
    { "sw_in", "SW_IN" }, // incoming shortwave flux, W m⁻²
    { "tair", "T_AIR" }, // Air temperature, ˚C
    { "taux", "IO_STRESS_X" }, // Ice-ocean stress x(east) direction, Pa
    { "tauy", "IO_STRESS_Y" }, // Ice-ocean stress y(north) direction, Pa
    { "tf", "TF" }, // Ocean freezing temperature, ˚C
    { "u", "ICE_U" }, // x(east)-ward ice velocity, m s⁻¹
    { "v", "ICE_V" }, // y(north)-ward ice velocity, m s⁻¹
    { "wind_speed", "WIND_SPEED" }, // Wind speed, m s⁻¹
    { "wind_u", "WIND_U" }, // wind velocity x component, m s⁻¹
    { "wind_v", "WIND_V" }, // wind velocity y component, m s⁻¹
