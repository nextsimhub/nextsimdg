/*!
 * @author  Tim Spain <timothy.spain@nersc.no>
 * @author  Einar Ólason <einar.olason@nersc.no>
 */

// External SharedArray names must be unique wrt to the external ProtectedArray names as well
{ "damage_upd", "DAMAGE" }, // Updated ice thickness, ice average, m
    { "hice_true", "H_ICE" }, // Updated ice thickness, ice average, m
    { "hice_dg", "H_ICE_DG" }, // Temporary DG hice tag
    { "cice_upd", "C_ICE" }, // Updated ice concentration
    { "cice_dg", "C_ICE_DG" }, // Temporary DG cice tag
    { "hsnow_true", "H_SNOW" }, // Updated snow depth, ice average, m
    { "tice_upd", "T_ICE" }, // Updated ice temperatures, ˚C
    { "qia", "Q_IA" }, // Ice to atmosphere heat flux W m⁻²
    { "qic", "Q_IC" }, // Ice conduction heat flux W m⁻²
    { "qio", "Q_IO" }, // Ice to ocean heat flux W m⁻²
    { "qow", "Q_OW" }, // Open water heat flux W m⁻²
    { "dqia_dt", "DQIA_DT" }, // Derivative of Qᵢₐ w.r.t. ice surface temperature  W m⁻² K⁻¹
    { "qpen_sw", "Q_PEN_SW" }, // Penetrating shortwave flux W m⁻²
    { "hsnow_melt", "HSNOW_MELT" }, // Thickness of snow that melted, m
    { "sublim", "SUBLIM" }, // Upward sublimation rate kg m⁻² s⁻¹
    { "delta_hice", "DELTA_HICE" }, // Change in sea ice thickness, m
    { "delta_cice", "DELTA_CICE" }, // Change in sea ice concentration
    { "new_ice", "NEW_ICE" }, // Volume of new ice formed [m]
    { "qsw_base", "Q_SW_BASE" }, // Short-wave flux through ice base W m⁻²
    { "qsw_ow", "Q_SW_OW" }, // Short-wave flux into ice free ocean W m⁻²
    { "taux_ow", "OW_STRESS_Y" }, // y(north)-ward open ocean stress, Pa
    { "tauy_ow", "OW_STRESS_X" }, // x(east)-ward open ocean stress, Pa
