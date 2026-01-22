/*!
 * @author  Tim Spain <timothy.spain@nersc.no>
 * @author Robert Jendersie <robert.jendersie@ovgu.de>
 */

#include "include/FiniteElementFluxes.hpp"

#include "include/Finalizer.hpp"
#include "include/FiniteElementSpecHum.hpp"
#include "include/IIceAlbedo.hpp"
#include "include/NextsimModule.hpp"
#include "include/constants.hpp"
#include "kokkos/include/KokkosTimer.hpp"

#include <memory>

namespace Nextsim {

double FiniteElementFluxes::dragOcean_q;
double FiniteElementFluxes::dragOcean_t;
double FiniteElementFluxes::dragIce_t;
double FiniteElementFluxes::m_oceanAlbedo;
double FiniteElementFluxes::m_I0;

static const double dragOcean_q_default = 1.5e-3;
static const double dragOcean_t_default = 0.83e-3;
static const double dragIce_t_default = 1.3e-3;
static const double oceanAlbedo_default = 0.07;
static const double i0_default = 0.17;

static const std::map<int, std::string> keyMap = {
    { FiniteElementFluxes::DRAGOCEANQ_KEY, "nextsim_thermo.drag_ocean_q" },
    { FiniteElementFluxes::DRAGOCEANT_KEY, "nextsim_thermo.drag_ocean_t" },
    { FiniteElementFluxes::DRAGICET_KEY, "nextsim_thermo.drag_ice_t" },
    { FiniteElementFluxes::OCEANALBEDO_KEY, "nextsim_thermo.albedoW" },
    { FiniteElementFluxes::I0_KEY, "nextsim_thermo.I_0" },
};

void FiniteElementFluxes::configure()
{
    Finalizer::registerUnique(Module::finalize<IIceAlbedo>);

    iIceAlbedoImpl = &Module::getImplementation<IIceAlbedo>();
    tryConfigure(iIceAlbedoImpl);

    dragOcean_q = Configured::getConfiguration(keyMap.at(DRAGOCEANQ_KEY), dragOcean_q_default);
    dragOcean_t = Configured::getConfiguration(keyMap.at(DRAGOCEANT_KEY), dragOcean_t_default);
    dragIce_t = Configured::getConfiguration(keyMap.at(DRAGICET_KEY), dragIce_t_default);
    m_oceanAlbedo = Configured::getConfiguration(keyMap.at(OCEANALBEDO_KEY), oceanAlbedo_default);
    m_I0 = Configured::getConfiguration(keyMap.at(I0_KEY), i0_default);
}

ConfigMap FiniteElementFluxes::getConfiguration() const
{
    return {
        { keyMap.at(DRAGOCEANQ_KEY), dragOcean_q },
        { keyMap.at(DRAGOCEANT_KEY), dragOcean_t },
        { keyMap.at(DRAGICET_KEY), dragIce_t },
        { keyMap.at(OCEANALBEDO_KEY), m_oceanAlbedo },
        { keyMap.at(I0_KEY), m_I0 },
    };
}

ModelState FiniteElementFluxes::getStateDiagnostic() const
{
    return { {
                 { "evap", evapAccessor.getHostRO() },
                 { "Q_lh_ow", Q_lh_owAccessor.getHostRO() },
                 { "Q_sh_ow", Q_sh_owAccessor.getHostRO() },
                 { "Q_lw_ow", Q_lw_owAccessor.getHostRO() },
                 { "Q_lh_ia", Q_lh_iaAccessor.getHostRO() },
                 { "Q_sh_ia", Q_sh_iaAccessor.getHostRO() },
                 { "Q_sw_ia", Q_sw_iaAccessor.getHostRO() },
                 { "Q_lw_ia", Q_lw_iaAccessor.getHostRO() },
                 { "rho_air", rho_airAccessor.getHostRO() },
                 { "cp_air", cp_airAccessor.getHostRO() },
                 { "sh_air", sh_airAccessor.getHostRO() },
                 { "sh_water", sh_waterAccessor.getHostRO() },
                 { "sh_ice", sh_iceAccessor.getHostRO() },
                 { "dshice_dT", dshice_dTAccessor.getHostRO() },
             },
        getConfiguration() };
}

void FiniteElementFluxes::setData(const ModelState::DataMap& ms)
{
    // Data arrays can now be set to the correct size
    Q_lh_owAccessor.getHostRW().resize();
    Q_sh_owAccessor.getHostRW().resize();
    Q_lw_owAccessor.getHostRW().resize();
    Q_lh_iaAccessor.getHostRW().resize();
    Q_sh_iaAccessor.getHostRW().resize();
    Q_sw_iaAccessor.getHostRW().resize();
    Q_lw_iaAccessor.getHostRW().resize();
    rho_airAccessor.getHostRW().resize();
    cp_airAccessor.getHostRW().resize();
    sh_airAccessor.getHostRW().resize();
    sh_waterAccessor.getHostRW().resize();
    sh_iceAccessor.getHostRW().resize();
    dshice_dTAccessor.getHostRW().resize();

    iIceAlbedoImpl->setData(ms);
}

FiniteElementFluxes::HelpMap& FiniteElementFluxes::getHelpText(HelpMap& map, bool getAll)
{
    map["FiniteElementFluxes"] = {
        { keyMap.at(DRAGOCEANQ_KEY), ConfigType::NUMERIC, { "0", "∞" },
            ConfigurationHelp::toString(dragOcean_q_default), "??",
            "Coefficient for evaporative mass flux calculation." },
        { keyMap.at(DRAGOCEANT_KEY), ConfigType::NUMERIC, { "0", "∞" },
            ConfigurationHelp::toString(dragOcean_t_default), "??",
            "Coefficient for sensible heat flux calculation." },
        { keyMap.at(DRAGICET_KEY), ConfigType::NUMERIC, { "0", "∞" },
            ConfigurationHelp::toString(dragIce_t_default), "??",
            "Ice drag coefficient for heat fluxes." },
        { keyMap.at(OCEANALBEDO_KEY), ConfigType::NUMERIC, { "0", "∞" },
            ConfigurationHelp::toString(oceanAlbedo_default), "",
            "Shortwave albedo of open ocean water." },
        { keyMap.at(I0_KEY), ConfigType::NUMERIC, { "0", "∞" },
            ConfigurationHelp::toString(i0_default), "", "Transmissivity of ice." },
    };
    return map;
}
FiniteElementFluxes::HelpMap& FiniteElementFluxes::getHelpRecursive(HelpMap& map, bool getAll)
{
    getHelpText(map, getAll);
    Module::getHelpRecursive<IIceAlbedo>(map, getAll);
    return map;
}

/*************************************************************/
// Drag coefficient from Gill(1982) / Smith (1980)
// Could be replaced by a  module ... but we'll probably never do that
KERNEL_IMPL_FUNCTION static double dragOcean_m(double windSpeed)
{
    return 1e-3 * Utils::max(1., Utils::min(2., 0.61 + 0.063 * windSpeed));
}

KERNEL_IMPL_FUNCTION static double latentHeatWater(double temperature)
{
    // Polynomial approximation expressed using Horner's scheme
    return Water::Lv0
        + temperature * (-2.36418e3 + temperature * (1.58927 + temperature * (-6.14342e-2)));
}

KERNEL_IMPL_FUNCTION static double latentHeatIce(double temperature)
{
    return Water::Lv0 + Water::Lf - 240. + temperature * (-290. + temperature * (-4.));
}

KERNEL_IMPL_FUNCTION static double stefanBoltzmannLaw(double temperatureC)
{
    return Ice::epsilon * PhysicalConstants::sigma * Utils::pow(kelvin(temperatureC), 4);
}

void FiniteElementFluxes::update(const TimestepTime& tst)
{
    static KokkosTimer<true> timer("FiniteElementFluxes::update");
    timer.start();
    updateAtmosphere(tst); // common atmospheric values
    updateOW(tst); // qow
    updateIce(tst); // qia & dqia/dT
    timer.stop();
}

/*************************************************************/
void FiniteElementFluxes::updateAtmosphere(const TimestepTime& tst)
{
    auto execSpace = DefaultExecutionSpace();
    auto& sh_water = sh_waterAccessor.getAutoRW(execSpace);
    auto& sh_ice = sh_iceAccessor.getAutoRW(execSpace);
    auto& cp_air = cp_airAccessor.getAutoRW(execSpace);
    auto& dshice_dT = dshice_dTAccessor.getAutoRW(execSpace);
    auto& rho_air = rho_airAccessor.getAutoRW(execSpace);
    auto& sh_air = sh_airAccessor.getAutoRW(execSpace);
    const auto& p_air = p_airAccessor.getAutoRO(execSpace);
    const auto& sss = sssAccessor.getAutoRO(execSpace);
    const auto& sst = sstAccessor.getAutoRO(execSpace);
    const auto& t_air = t_airAccessor.getAutoRO(execSpace);
    const auto& t_dew2 = t_dew2Accessor.getAutoRO(execSpace);
    const auto& tsurf = tsurfAccessor.getAutoRO(execSpace);

    const FiniteElementSpecHum specHumWater = FiniteElementSpecHum::water();
    const FiniteElementSpecHum specHumIce = FiniteElementSpecHum::ice();

    static KokkosTimer<true> timer("FiniteElementFluxes::updateAtmosphere");
    timer.start();
    overElementsAuto(OVER_ELEMENTS_LAMBDA(const ElementIndex i) {
        // Specific humidity of...
        // ...the air
        sh_air[i] = specHumWater(t_dew2[i], p_air[i]);
        // ...over the open ocean
        sh_water[i] = specHumWater(sst[i], p_air[i], sss[i]);
        // ...over the ice
        const auto [shIce, dshIceDT] = specHumIce.valueAndDerivative(tsurf[i], p_air[i]);
        sh_ice[i] = shIce;
        dshice_dT[i] = dshIceDT;
        // Density of the wet air
        double Ra_wet = Air::Ra / (1 - sh_air[i] * (1 - Vapour::Ra / Air::Ra));
        rho_air[i] = p_air[i] / (Ra_wet * kelvin(t_air[i]));
        // Heat capacity of the wet air
        cp_air[i] = Air::cp + sh_air[i] * Vapour::cp;
    });
    timer.stop();
}

void FiniteElementFluxes::updateOW(const TimestepTime& tst)
{
    auto execSpace = DefaultExecutionSpace();
    auto& Q_lh_ow = Q_lh_owAccessor.getAutoRW(execSpace);
    auto& evap = evapAccessor.getAutoRW(execSpace);
    auto& qow = qowAccessor.getAutoRW(execSpace);
    auto& Q_sw_ow = Q_sw_owAccessor.getAutoRW(execSpace);
    auto& Q_sh_ow = Q_sh_owAccessor.getAutoRW(execSpace);
    auto& tau_x_ow = tau_x_owAccessor.getAutoRW(execSpace);
    auto& Q_lw_ow = Q_lw_owAccessor.getAutoRW(execSpace);
    auto& tau_y_ow = tau_y_owAccessor.getAutoRW(execSpace);
    const auto& cp_air = cp_airAccessor.getAutoRO(execSpace);
    const auto& lw_in = lw_inAccessor.getAutoRO(execSpace);
    const auto& rho_air = rho_airAccessor.getAutoRO(execSpace);
    const auto& sh_air = sh_airAccessor.getAutoRO(execSpace);
    const auto& sh_water = sh_waterAccessor.getAutoRO(execSpace);
    const auto& sst = sstAccessor.getAutoRO(execSpace);
    const auto& sw_in = sw_inAccessor.getAutoRO(execSpace);
    const auto& t_air = t_airAccessor.getAutoRO(execSpace);
    const auto& u_air = u_airAccessor.getAutoRO(execSpace);
    const auto& v_air = v_airAccessor.getAutoRO(execSpace);
    const auto& windSpeed = windSpeedAccessor.getAutoRO(execSpace);

    // static members can not be captured directly
    const double dragOcean_q = FiniteElementFluxes::dragOcean_q;
    const double dragOcean_t = FiniteElementFluxes::dragOcean_t;
    const double m_oceanAlbedo = FiniteElementFluxes::m_oceanAlbedo;

    overElementsAuto(OVER_ELEMENTS_LAMBDA(const ElementIndex i) {
        // Mass flux from open water (evaporation)
        evap[i] = dragOcean_q * rho_air[i] * windSpeed[i] * (sh_water[i] - sh_air[i]);
        // Momentum flux from open water (drag pressure)
        /* Drag the ocean experiences from the wind - still only used in the coupled case */
        const double oceanDrag = dragOcean_m(windSpeed[i]);
        tau_x_ow[i] = rho_air[i] * oceanDrag * u_air[i] * windSpeed[i];
        tau_y_ow[i] = rho_air[i] * oceanDrag * v_air[i] * windSpeed[i];

        // Heat flux open water
        //   Latent heat from evaporation (and condensation)
        Q_lh_ow[i] = evap[i] * latentHeatWater(sst[i]);
        //   Sensible heat
        Q_sh_ow[i] = dragOcean_t * rho_air[i] * cp_air[i] * windSpeed[i] * (sst[i] - t_air[i]);
        //   Shortwave flux
        Q_sw_ow[i] = -sw_in[i] * (1 - m_oceanAlbedo);
        // Longwave flux
        Q_lw_ow[i] = stefanBoltzmannLaw(sst[i]) - lw_in[i];
        // Total open water flux
        qow[i] = Q_lh_ow[i] + Q_sh_ow[i] + Q_sw_ow[i] + Q_lw_ow[i];
    });
}

void FiniteElementFluxes::updateIce(const TimestepTime& tst)
{
    auto execSpace = DefaultExecutionSpace();
    auto& qia = qiaAccessor.getAutoRW(execSpace);
    auto& subl = sublAccessor.getAutoRW(execSpace);
    auto& Q_lh_ia = Q_lh_iaAccessor.getAutoRW(execSpace);
    auto& Q_lw_ia = Q_lw_iaAccessor.getAutoRW(execSpace);
    auto& penSW = penSWAccessor.getAutoRW(execSpace);
    auto& dqia_dt = dqia_dtAccessor.getAutoRW(execSpace);
    auto& Q_sw_base = Q_sw_baseAccessor.getAutoRW(execSpace);
    auto& Q_sw_ia = Q_sw_iaAccessor.getAutoRW(execSpace);
    auto& Q_sh_ia = Q_sh_iaAccessor.getAutoRW(execSpace);
    const auto& cice = ciceAccessor.getAutoRO(execSpace);
    const auto& cp_air = cp_airAccessor.getAutoRO(execSpace);
    const auto& dshice_dT = dshice_dTAccessor.getAutoRO(execSpace);
    const auto& hsnow = hsnowAccessor.getAutoRO(execSpace);
    const auto& lw_in = lw_inAccessor.getAutoRO(execSpace);
    const auto& rho_air = rho_airAccessor.getAutoRO(execSpace);
    const auto& sh_air = sh_airAccessor.getAutoRO(execSpace);
    const auto& sh_ice = sh_iceAccessor.getAutoRO(execSpace);
    const auto& sw_in = sw_inAccessor.getAutoRO(execSpace);
    const auto& t_air = t_airAccessor.getAutoRO(execSpace);
    const auto& tsurf = tsurfAccessor.getAutoRO(execSpace);
    const auto& windSpeed = windSpeedAccessor.getAutoRO(execSpace);
    const auto& iceAlbedo = iceAlbedoAccessor.getAutoRO(execSpace);
    const auto& icePenSW = icePenSWAccessor.getAutoRO(execSpace);

    // static members can not be captured directly
    const double dragIce_t = FiniteElementFluxes::dragIce_t;

    iIceAlbedoImpl->setTime(tst.start);
    iIceAlbedoImpl->update(tst);
    overElementsAuto(OVER_ELEMENTS_LAMBDA(const ElementIndex i) {
        // Mass flux ice
        subl[i] = dragIce_t * rho_air[i] * windSpeed[i] * (sh_ice[i] - sh_air[i]);

        // Momentum flux is dealt with by the ice dynamics

        // Heat flux ice-atmosphere
        // Latent heat from sublimation
        Q_lh_ia[i] = subl[i] * latentHeatIce(tsurf[i]);
        double dmdot_dT = dragIce_t * rho_air[i] * windSpeed[i] * dshice_dT[i];
        double dQlh_dT = latentHeatIce(tsurf[i]) * dmdot_dT;

        // Sensible heat flux
        Q_sh_ia[i] = dragIce_t * rho_air[i] * cp_air[i] * windSpeed[i] * (tsurf[i] - t_air[i]);
        double dQsh_dT = dragIce_t * rho_air[i] * cp_air[i] * windSpeed[i];

        // Shortwave flux
        const double albedoValue = iceAlbedo[i];
        const double i0 = icePenSW[i];
        Q_sw_ia[i] = -sw_in[i] * (1. - albedoValue) * (1. - i0);
        const double extinction = 0.; // TODO: Replace with de Beer's law or a module
        penSW[i] = sw_in[i] * (1. - albedoValue) * i0 * (1. - extinction);
        Q_sw_base[i] = sw_in[i] * (1. - albedoValue) * i0 * extinction;

        // Longwave flux
        Q_lw_ia[i] = stefanBoltzmannLaw(tsurf[i]) - lw_in[i];
        double dQlw_dT = 4 / kelvin(tsurf[i]) * stefanBoltzmannLaw(tsurf[i]);

        // Total flux
        qia[i] = Q_lh_ia[i] + Q_sh_ia[i] + Q_sw_ia[i] + Q_lw_ia[i];

        // Total temperature dependence of flux
        dqia_dt[i] = dQlh_dT + dQsh_dT + dQlw_dT;
    });
}

} /* namespace Nextsim */
