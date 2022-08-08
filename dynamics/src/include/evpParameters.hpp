/*!
 * @file evpParameters.hpp
 * @date August 5 2022
 * @author Thomas Richter <thomas.richter@ovgu.de>
 */

#ifndef __EVPPARAMETERS_HPP
#undef __EVPPARAMETERS_HPP

namespace Nextsim {

class EVPParameters {
public:
    constexpr double rho_ice = 900.0; //!< Sea ice density
    constexpr double rho_atm = 1.3; //!< Air density
    constexpr double rho_ocean = 1026.0; //!< Ocean density

    constexpr double C_atm = 1.2e-3; //!< Air drag coefficient
    constexpr double C_ocean = 5.5e-3; //!< Ocean drag coefficient

    constexpr double F_atm = C_atm * rho_atm; //!< effective factor for atm-forcing
    constexpr double F_ocean = C_ocean * rho_ocean; //!< effective factor for ocean-forcing

    constexpr double Pstar = 27500.0; //!< Ice strength
    constexpr double fc = 1.46e-4; //!< Coriolis

    constexpr double DeltaMin = 2.e-9; //!< Viscous regime
};

} /* namespace Nextsim */

#endif /* __MEVP_HPP */
