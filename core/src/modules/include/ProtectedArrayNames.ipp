/*!
 * @file ProtectedArrayNames.ipp
 *
 * @date 1 Jul 2024
 * @author Tim Spain <timothy.spain@nersc.no>
 * @author Einar Ólason <einar.olason@nersc.no>
 */

// External ProtectedArray names must be unique wrt to the external SharedArray names as well
{ "tair", "T_AIR" }, // Air temperature, ˚C
{ "dew2m", "DEW_2M" }, // Dew point at 2 m, ˚C
{ "pair", "P_AIR" }, // sea level air pressure, Pa
{ "mixrat", "MIXRAT" }, // water vapour mass mixing ratio
{ "sw_in", "SW_IN" }, // incoming shortwave flux, W m⁻²
{ "lw_in", "LW_IN" }, // incoming longwave flux, W m⁻²
{ "mld", "MLD" }, // mixed layer depth, m
{ "snowfall", "SNOWFALL" }, // snow fall, kg m⁻² s⁻¹
{ "sss", "SSS" }, // sea surface salinity, PSU
{ "sst", "SST" }, // sea surface temperature ˚C
{ "sst_ext", "EXT_SST" }, // External sea surface temperature ˚C
{ "sss_ext", "EXT_SSS" }, // External sea surface salinity PSU
{ "eminusp", "E-P" }, // E-P atmospheric freshwater flux, kg s⁻¹ m⁻²
{ "mlcp", "CPML" }, // Mixed layer bulk heat capacity J K⁻¹ m⁻²
{ "tf", "TF" }, // Ocean freezing temperature, ˚C
{ "wind_speed", "WIND_SPEED" }, // Wind speed, m s⁻¹
{ "wind_u", "WIND_U" }, // wind velocity x component, m s⁻¹
{ "wind_v", "WIND_V" }, // wind velocity y component, m s⁻¹
{ "hice_true_pro", "HTRUE_ICE" }, // Ice thickness, ice average, m
{ "hsnow_true_pro", "HTRUE_SNOW" }, // Snow thickness, ice average, m
{ "ocean_u", "OCEAN_U" }, // x(east)-ward ocean current, m s⁻¹
{ "ocean_v", "OCEAN_V" }, // y(north)-ward ocean current, m s⁻¹
{ "u", "ICE_U" }, // x(east)-ward ice velocity, m s⁻¹
{ "v", "ICE_V" }, // y(north)-ward ice velocity, m s⁻¹
{ "sst_slab", "SLAB_SST" }, // Slab ocean surface temperature ˚C
{ "sss_slab", "SLAB_SSS" }, // Slab ocean surface salinity PSU
{ "qdw", "SLAB_QDW" }, // Slab ocean temperature nudging heat flux, W m⁻²
{ "fdw", "SLAB_FDW" }, // Slab ocean salinity nudging water flux, kg s⁻¹ m⁻²
{ "e11", "E11" }, // The 1 1 component of strain rate, s⁻¹
{ "e22", "E22" }, // The 2 2 component of strain rate, s⁻¹
{ "e12", "E12" }, // The 1 2 component of strain rate, s⁻¹
{ "s11", "S11" }, // The 1 1 component of stress, Pa
{ "s22", "S22" }, // The 2 2 component of stress, Pa
{ "s12", "S12" }, // The 1 2 component of stress, Pa

