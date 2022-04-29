/*!
 * @file FiniteElementFluxes.hpp
 *
 * @date Apr 29, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef FINITEELEMENTFLUXES_HPP
#define FINITEELEMENTFLUXES_HPP

#include "include/IFluxCalculation.hpp"

#include "include/ModelArrayRef.hpp"

namespace Nextsim {

class FiniteElementFluxes : public IFluxCalculation {
public:
    FiniteElementFluxes()
        : IFluxCalculation()
    {
    }

    ModelState getState() const override;
    ModelState getState(const OutputLevel&) const override;

    std::string getName() const override { return "FiniteElementFluxes"; }

    void update(const TimestepTime& tst) override;

private:
    // Owned diagnostic fields
    HField evap; // Open water evaporative mass flux [kg  m⁻²]
    HField Q_lh_ow; // Open water latent heat flux [W m⁻²]
    HField Q_sh_ow; // Open water sensible heat flux [W m⁻²]
    HField Q_sw_ow; // Open water incident shortwave radiative flux [W m⁻²]
    HField Q_lw_ow; // Open water net longwave radiative flux [W m⁻²]
    HField Q_lh_i; // Ice latent heat flux [W m⁻²]
    HField Q_sh_i; // Ice sensible heat flux [W m⁻²]
    HField Q_sw_i; // Ice incident shortwave radiative flux [W m⁻²]
    HField Q_lw_i; // Ice net longwave radiative flux [W m⁻²]
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
    ModelArrayRef<ProtectedArray::SW_IN> sw_in;
    ModelArrayRef<ProtectedArray::LW_IN> lw_in;


};

} /* namespace Nextsim */

#endif /* FINITEELEMENTFLUXES_HPP */
