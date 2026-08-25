/*!
 * @author  Tim Spain <timothy.spain@nersc.no>
 */

#ifndef SRC_INCLUDE_CONSTANTS_HPP
#define SRC_INCLUDE_CONSTANTS_HPP

#include "include/FloatType.hpp"
#include "include/KernelAlternatives.hpp"

namespace Nextsim {

//! General physical constants of the Earth and universe
namespace PhysicalConstants {

    /*!
     *  Standard acceleration due to gravity at the Earth's poles. [m s⁻²]
     *
     *  WGS 84 Ellipsoidal gravity formula evaluated at 90˚ latitude
     */
    const FloatType g = 9.8321849378;

    //! Stefan-Boltzmann constant. [W m⁻² K⁻⁴]
    const FloatType sigma = 5.670374419e-8;

    //! Von Karman constant. [1]
    const FloatType vonKarman = 0.4;

    //! Rotation rate of the Earth [rad s⁻¹]
    const FloatType omega = 7.2921158e-5;

    //! Triple point temperature of pure water [K]
    const FloatType Tt = 273.16;

    //! Ratio of circumference to radius
    const FloatType tau = 6.28318530717958647652;

    // degrees to radians as a hex float
    const FloatType deg2rad = 0x1.1df46a2529d39p-6;

    // radians to degrees
    const FloatType rad2deg = 0x1.ca5dc1a63c1f8p+5;
}

//! Properties of water ice around 0˚C and 101.3 kPa
namespace Ice {

    //! Specific heat capacity at constant pressure of water ice [J kg⁻¹ K⁻¹]
    const FloatType cp = 2100.;

    //! Thermal emissivity of smooth ice [0..1]
    const FloatType epsilon = 0.996;

    //! Heat conductivity of ice [W m⁻¹ K⁻¹]
    const FloatType kappa = 2.0334;

    //! Latent heat of fusion of ice/water [J kg⁻¹]
    const FloatType Lf = 333.55e3;

    /*!
     * Density of ice. [kg m⁻³]
     *
     * Taken to be the same value as used in NEMO-LIM.
     */
    const FloatType rho = 917;

    /*!
     * Density of snow. [kg m⁻³]
     *
     * Taken to be the same value as used in NEMO-LIM.
     */
    const FloatType rhoSnow = 330.;

    //! Salinity of sea ice. [g kg⁻¹]
    const FloatType s = 5;

    //! Melting point of pure ice [K]
    const FloatType Tm = 273.15;
}

//! Properties of dry air around 0˚C and 101.3 kPa
namespace Air {

    //! Specific heat capacity at constant pressure of dry air [J kg⁻¹ K⁻¹]
    const FloatType cp = 1004.64;

    //! Specific gas constant for dry air [J kg⁻¹ K⁻¹]
    const FloatType Ra = 287.058;

    //! Density of dry air at IUPAC STP [kg m⁻³]
    const FloatType rho = 1.2754;
}

//! Properties of water vapour
namespace Vapour {

    //! Specific heat capacity at constant pressure of water vapour [J kg⁻¹ K⁻¹]
    const FloatType cp = 1860.;

    //! Latent heat of vaporization at 0˚C [J kg⁻¹]
    const FloatType Lv0 = 2500.79e3;

    //! Specific gas constant for water vapour [J kg⁻¹ K⁻¹]
    const FloatType Ra = 461.5;
}

//! Properties of liquid water
namespace Water {

    //! Specific heat capacity at constant pressure of water [J kg⁻¹ K⁻¹]
    const FloatType cp = 4186.84;

    //! Latent heat of fusion of water/ice [J kg⁻¹]
    const FloatType Lf = Ice::Lf;

    //! Latent heat of vaporization at 0˚C [J kg⁻¹]
    const FloatType Lv0 = Vapour::Lv0;

    //! Proportionality constant between salinity in psu and freezing point depression [K psu⁻¹]
    const FloatType mu = 0.055;

    //! Density of fresh water at 4˚C. [kg m⁻³]
    const FloatType rho = 1000;

    //! Typical density of ocean water. [kg m⁻³]
    const FloatType rhoOcean = 1025;

    //! Freezing point of pure water. [K]
    const FloatType Tf = Ice::Tm;

    //! Freezing point of typical ocean water. [˚C]
    const FloatType TfOcean = -1.8;
}

//! Convert a temperature from ˚C to K
KERNEL_IMPL_FUNCTION inline FloatType kelvin(FloatType celsius) { return celsius + Water::Tf; }

//! Convert a temperature from K to ˚C
KERNEL_IMPL_FUNCTION inline FloatType celsius(FloatType kelvin) { return kelvin - Water::Tf; }

//! Convert an angle from radians to degrees
KERNEL_IMPL_FUNCTION inline FloatType degrees(FloatType radians)
{
    return radians * PhysicalConstants::rad2deg;
}

//! Convert an angle from degrees to radians
KERNEL_IMPL_FUNCTION inline FloatType radians(FloatType degrees)
{
    return degrees * PhysicalConstants::deg2rad;
}

//! Convert a pressure from Pa to mbar
KERNEL_IMPL_FUNCTION inline FloatType mbar(FloatType pascals) { return pascals / 100; }

//! Convert a pressure from mbar to Pa
KERNEL_IMPL_FUNCTION inline FloatType pascals(FloatType mbar) { return mbar * 100; }

} /* namespace Nextsim */

#endif /* SRC_INCLUDE_CONSTANTS_HPP */
