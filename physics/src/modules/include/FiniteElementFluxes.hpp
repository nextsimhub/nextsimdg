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
#include "include/IIceFluxes.hpp"
#include "include/IOWFluxes.hpp"

#include "include/ModelArrayRef.hpp"

namespace Nextsim {

class FiniteElementFluxes : public IIceFluxes, public IOWFluxes, public Configured<FiniteElementFluxes> {
public:
    FiniteElementFluxes()
        : IIceFluxes()
        , IOWFluxes()
    {
    }
    ~FiniteElementFluxes() = default;

    ModelState getState() const override;
    ModelState getState(const OutputLevel&) const override;

    std::string getName() const override { return "FiniteElementFluxes"; }

    void updateOW(const TimestepTime& tst) override
    {
        overElements(std::bind(&FiniteElementFluxes::calculateOW, this, std::placeholders::_1,
                         std::placeholders::_2),
            tst);
    }

    void updateIce(const TimestepTime& tst) override
    {
        overElements(std::bind(&FiniteElementFluxes::calculateIce, this, std::placeholders::_1,
                         std::placeholders::_2),
            tst);
    }

    void updateAtmosphere(const TimestepTime& tst);

    void configure() override;
    enum {
        DRAGOCEANQ_KEY,
        DRAGOCEANT_KEY,
        DRAGICET_KEY,
        OCEANALBEDO_KEY,
        I0_KEY,
    };

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
    // Input fields
    ModelArrayRef<ProtectedArray::SST> sst;
    ModelArrayRef<ProtectedArray::T_AIR> t_air;
    ModelArrayRef<ProtectedArray::P_AIR> p_air;
    ModelArrayRef<ProtectedArray::DENSITY_AIR> rho_air;
    ModelArrayRef<ProtectedArray::WIND_SPEED> v_air;
    ModelArrayRef<ProtectedArray::CP_AIR> cp_air;
    ModelArrayRef<ProtectedArray::SPEC_HUM_AIR> sh_air;
    ModelArrayRef<ProtectedArray::SPEC_HUM_WATER> sh_water;
    ModelArrayRef<ProtectedArray::SPEC_HUM_ICE> sh_ice;
    ModelArrayRef<ProtectedArray::H_SNOW> h_snow; // cell-averaged value
    ModelArrayRef<ProtectedArray::C_ICE> cice;
    ModelArrayRef<ProtectedArray::T_ICE> tice;
    ModelArrayRef<ProtectedArray::SW_IN> sw_in;
    ModelArrayRef<ProtectedArray::LW_IN> lw_in;

    void calculateOW(size_t i, const TimestepTime& tst);
    void calculateIce(size_t i, const TimestepTime& tst);

    static double dragOcean_q;
    static double dragOcean_m(double windSpeed);
    static double dragOcean_t;
    static double dragIce_t;

    static double m_oceanAlbedo;

    static double m_I0;

    static double latentHeatWater(double temperature);
    static double latentHeatIce(double temperature);

    std::unique_ptr<IIceAlbedo> iIceAlbedoImpl;
};

class FiniteElementFluxCalc : public IFluxCalculation, public Configured<FiniteElementFluxCalc> {
public:
    FiniteElementFluxCalc()
    : IFluxCalculation()
    , fef(nullptr)
    , iIceFluxesImpl(nullptr)
    {
    }

    void setData(const ModelState&) override { }

    ModelState getState() const override;
    ModelState getState(const OutputLevel&) const override;

    std::string getName() const override { return "FiniteElementFluxCalc"; }

    void update(const TimestepTime&) override;
    enum {
        OW_FLUX_KEY,
    };
    void configure() override;

private:
    std::unique_ptr<IOWFluxes> iOWFluxesImpl;
    IIceFluxes* iIceFluxesImpl;
    FiniteElementFluxes* fef;
};

} /* namespace Nextsim */

#endif /* FINITEELEMENTFLUXES_HPP */
