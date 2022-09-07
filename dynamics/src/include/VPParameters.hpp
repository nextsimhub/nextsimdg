/*!
 * @file VPParameters.hpp
 * @date August 5 2022
 * @author Thomas Richter <thomas.richter@ovgu.de>
 */

#ifndef __VPPARAMETERS_HPP
#define __VPPARAMETERS_HPP

namespace Nextsim {

class VPParameters {
public:
    double rho_ice; //!< Sea ice density
    double rho_atm; //!< Air density
    double rho_ocean; //!< Ocean density

    double C_atm; //!< Air drag coefficient
    double C_ocean; //!< Ocean drag coefficient

    double F_atm; //!< effective factor for atm-forcing
    double F_ocean; //!< effective factor for ocean-forcing

    double Pstar; //!< Ice strength
    double fc; //!< Coriolis

    double DeltaMin; //!< Viscous regime

    VPParameters()
    {
        rho_ice = 900.0; //!< Sea ice density
        rho_atm = 1.3; //!< Air density
        rho_ocean = 1026.0; //!< Ocean density

        C_atm = 1.2e-3; //!< Air drag coefficient
        C_ocean = 5.5e-3; //!< Ocean drag coefficient

        F_atm = C_atm * rho_atm; //!< effective factor for atm-forcing
        F_ocean = C_ocean * rho_ocean; //!< effective factor for ocean-forcing

        Pstar = 27500.0; //!< Ice strength
        fc = 1.46e-4; //!< Coriolis

        DeltaMin = 2.e-9; //!< Viscous regime
    }
};

} /* namespace Nextsim */

#endif /* __VP_HPP */
