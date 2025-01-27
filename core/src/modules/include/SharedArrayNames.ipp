/*!
 * @file SharedArrayNames.ipp
 *
 * @date 1 Jul 2024
 * @author Tim Spain <timothy.spain@nersc.no>
 * @author Einar Ólason <einar.olason@nersc.no>
 */

// External SharedArray names must be unique wrt to the external ProtectedArray names as well
{ "cice", "C_ICE" }, // Updated ice concentration
    { "damage", "DAMAGE" }, // Updated ice thickness, ice average, m
    { "delta_cice", "DELTA_CICE" }, // Change in sea ice concentration
    { "delta_hice", "DELTA_HICE" }, // Change in sea ice thickness, m
    { "dqia_dt", "DQIA_DT" }, // Derivative of Qᵢₐ w.r.t. ice surface temperature  W m⁻² K⁻¹
    { "hice", "H_ICE" }, // Updated ice thickness, ice average, m
    { "hsnow", "H_SNOW" }, // Updated snow depth, ice average, m
    { "hsnow_melt", "HSNOW_MELT" }, // Thickness of snow that melted, m
    { "new_ice", "NEW_ICE" }, // Volume of new ice formed [m]
    { "qia", "Q_IA" }, // Ice to atmosphere heat flux W m⁻²
    { "qic", "Q_IC" }, // Ice conduction heat flux W m⁻²
    { "qio", "Q_IO" }, // Ice to ocean heat flux W m⁻²
    { "qow", "Q_OW" }, // Open water heat flux W m⁻²
    { "qpen_sw", "Q_PEN_SW" }, // Penetrating shortwave flux W m⁻²
    { "sublim", "SUBLIM" }, // Upward sublimation rate kg m⁻² s⁻¹
    { "tice", "T_ICE" }, // Updated ice temperatures, ˚C
