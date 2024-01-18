/*!
 * @file CommonDynamicsParameters.hpp
 *
 * @date Jan 18, 2024
 * @author Tim Spain <timothy.spain@nersc.no>
 * @author Thomas Richter <thomas.richter@ovgu.de>
 *
 */

#ifndef COMMONDYNAMICSPARAMETERS_HPP
#define COMMONDYNAMICSPARAMETERS_HPP

namespace Nextsim {
class CommonDynamicsParameters {
public:
    double rho_ice; //!< Sea ice density
    double rho_atm; //!< Air density
    double rho_ocean; //!< Ocean density

    double C_atm; //!< Air drag coefficient
    double C_ocean; //!< Ocean drag coefficient

    double F_atm; //!< effective factor for atm-forcing
    double F_ocean; //!< effective factor for ocean-forcing

    double fc; //!< Coriolis

    CommonDynamicsParameters()
    {
        rho_ice = 900.0; //!< Sea ice density
        rho_atm = 1.3; //!< Air density
        rho_ocean = 1026.0; //!< Ocean density

        C_atm = 1.2e-3; //!< Air drag coefficient
        C_ocean = 5.5e-3; //!< Ocean drag coefficient

        F_atm = C_atm * rho_atm; //!< effective factor for atm-forcing
        F_ocean = C_ocean * rho_ocean; //!< effective factor for ocean-forcing

        fc = 1.46e-4; //!< Coriolis
    }
};
}

#endif /* COMMONDYNAMICSPARAMETERS_HPP */
