/*!
 * @file FiniteElementFluxes.hpp
 *
 * @date Apr 29, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef FINITEELEMENTFLUXES_HPP
#define FINITEELEMENTFLUXES_HPP

#include "include/Configured.hpp"
#include "include/IFluxCalculation.hpp"
#include "include/IIceAlbedo.hpp"
#include "include/IIceOceanHeatFlux.hpp"
#include "include/ModelArrayRef.hpp"
#include "include/ModelComponent.hpp"

namespace Nextsim {

//! A class that implements ice fluxes and open water fluxes according finiteelement.cpp.
class FiniteElementFluxes : public IFluxCalculation, public Configured<FiniteElementFluxes> {
public:
    FiniteElementFluxes()
        : iIceAlbedoImpl(nullptr)
        , evap(ModelArray::Type::H)
        , Q_lh_ow(ModelArray::Type::H)
        , Q_sh_ow(ModelArray::Type::H)
        , Q_sw_ow(ModelArray::Type::H)
        , Q_lw_ow(ModelArray::Type::H)
        , Q_lh_ia(ModelArray::Type::H)
        , Q_sh_ia(ModelArray::Type::H)
        , Q_sw_ia(ModelArray::Type::H)
        , Q_lw_ia(ModelArray::Type::H)
        , rho_air(ModelArray::Type::H)
        , cp_air(ModelArray::Type::H)
        , sh_air(ModelArray::Type::H)
        , sh_water(ModelArray::Type::H)
        , sh_ice(ModelArray::Type::H)
        , dshice_dT(ModelArray::Type::H)
        , sst(getStore())
        , sss(getStore())
        , t_air(getStore())
        , t_dew2(getStore())
        , p_air(getStore())
        , v_air(getStore())
        , h_snow(getStore())
        , h_snow_true(getStore())
        , cice(getStore())
        , tice(getStore())
        , sw_in(getStore())
        , lw_in(getStore())
    {
    }
    ~FiniteElementFluxes() = default;

    enum {
        DRAGOCEANQ_KEY,
        DRAGOCEANT_KEY,
        DRAGICET_KEY,
        OCEANALBEDO_KEY,
        I0_KEY,
    };
    void configure() override;

    void setData(const ModelState::DataMap&) override;

    ModelState getState() const override;
    ModelState getState(const OutputLevel&) const override;
    ModelState getStateRecursive(const OutputSpec& os) const override;

    static HelpMap& getHelpText(HelpMap& map, bool getAll);
    static HelpMap& getHelpRecursive(HelpMap& map, bool getAll);

    std::string getName() const override { return "FiniteElementFluxes"; }

    void update(const TimestepTime& tst) override;

    //! Updates the fluxes over open water.
    void updateOW(const TimestepTime& tst);

    //! Updates the fluxes over ice.
    void updateIce(const TimestepTime& tst);

    //! Updates the atmospheric fluxes.
    void updateAtmosphere(const TimestepTime& tst);

private:
    // Owned diagnostic fields
    HField evap; // Open water evaporative mass flux [kg  m⁻²]
    HField Q_lh_ow; // Open water latent heat flux [W m⁻²]
    HField Q_sh_ow; // Open water sensible heat flux [W m⁻²]
    HField Q_sw_ow; // Open water incident shortwave radiative flux [W m⁻²]
    HField Q_lw_ow; // Open water net longwave radiative flux [W m⁻²]
    HField Q_lh_ia; // Ice latent heat flux [W m⁻²]
    HField Q_sh_ia; // Ice sensible heat flux [W m⁻²]
    HField Q_sw_ia; // Ice incident shortwave radiative flux [W m⁻²]
    HField Q_lw_ia; // Ice net longwave radiative flux [W m⁻²]
    // Derived air properties
    HField rho_air;
    HField cp_air; // Specific heat capacity of the wet air [J kg⁻¹ K⁻¹]
    // Specific humidity and T derivative
    HField sh_air;
    HField sh_water;
    HField sh_ice;
    HField dshice_dT;
    // Input fields
    ModelArrayRef<Protected::SST> sst;
    ModelArrayRef<Protected::SSS> sss;
    ModelArrayRef<Protected::T_AIR> t_air;
    ModelArrayRef<Protected::DEW_2M> t_dew2;
    ModelArrayRef<Protected::P_AIR> p_air;
    ModelArrayRef<Protected::WIND_SPEED> v_air;
    ModelArrayRef<Protected::H_SNOW> h_snow; // cell-averaged value
    ModelArrayRef<Protected::HTRUE_SNOW> h_snow_true; // cell-averaged value
    ModelArrayRef<Protected::C_ICE> cice;
    ModelArrayRef<Protected::T_ICE> tice;
    ModelArrayRef<Protected::SW_IN> sw_in;
    ModelArrayRef<Protected::LW_IN> lw_in;

    void calculateOW(size_t i, const TimestepTime& tst);
    void calculateIce(size_t i, const TimestepTime& tst);
    void calculateAtmos(size_t i, const TimestepTime& tst);

    static double dragOcean_q;
    static double dragOcean_m(double windSpeed);
    static double dragOcean_t;
    static double dragIce_t;

    static double m_oceanAlbedo;

    static double m_I0;

    static double latentHeatWater(double temperature);
    static double latentHeatIce(double temperature);

    IIceAlbedo* iIceAlbedoImpl;
};

} /* namespace Nextsim */

#endif /* FINITEELEMENTFLUXES_HPP */
