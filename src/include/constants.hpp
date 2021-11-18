/*!
 * @file constants.hpp
 * @date Sep 14, 2021
 * @author Tim Spain
 */

#ifndef SRC_INCLUDE_CONSTANTS_HPP
#define SRC_INCLUDE_CONSTANTS_HPP

//! General physical constants of the Earth and universe
namespace PhysicalConstants {

/*!
 *  Standard acceleration due to gravity at the Earth's poles. [m s⁻²]
 *
 *  WGS 84 Ellipsoidal gravity formula evaluated at 90˚ latitude
 */
const double g = 9.8321849378;

//! Stefan-Boltzmann constant. [W m⁻² K⁻⁴]
const double sigma = 5.670374419e-8;

//! Von Karman constant. [1]
const double vonKarman = 0.4;

//! Rotation rate of the Earth [rad s⁻¹]
const double omega = 7.2921158e-5;

//! Triple point temperature of pure water [K]
const double Tt = 273.16;

//! Ratio of circumference to radius
const double tau = 6.28318530717958647652;
}

//! Properties of water ice around 0˚C and 101.3 kPa
namespace Ice {

//! Specific heat capacity at constant pressure of water ice [J kg⁻¹ K⁻¹]
const double cp = 2100.;

//! Thermal emissivity of smooth ice [0..1]
const double epsilon = 0.996;

//! Heat conductivity of ice [W m⁻¹ K⁻¹]
const double kappa = 2.0334;

//! Latent heat of fusion of ice/water [J kg⁻¹]
const double Lf = 333.55e3;

/*!
 * Density of ice. [kg m⁻³]
 *
 * Taken to be the same value as used in NEMO-LIM.
 */
const double rho = 917;

/*!
 * Density of snow. [kg m⁻³]
 *
 * Taken to be the same value as used in NEMO-LIM.
 */
const double rhoSnow = 330.;

//! Salinity of sea ice. [g kg⁻¹]
const double s = 5;

//! Melting point of pure ice [K]
const double Tm = 273.15;
}

//! Properties of dry air around 0˚C and 101.3 kPa
namespace Air {

//! Specific heat capacity at constant pressure of dry air [J kg⁻¹ K⁻¹]
const double cp = 1004.64;

//! Specific gas constant for dry air [J kg⁻¹ K⁻¹]
const double Ra = 287.058;

//! Density of dry air at IUPAC STP [kg m⁻³]
const double rho = 1.2754;
}

//! Properties of water vapour
namespace Vapour {

//! Specific heat capacity at constant pressure of water vapour [J kg⁻¹ K⁻¹]
const double cp = 1860.;

//! Latent heat of vaporization at 0˚C [J kg⁻¹]
const double Lv0 = 2500.79e3;

//! Specific gas constant for water vapour [J kg⁻¹ K⁻¹]
const double Ra = 461.5;
}

//! Properties of liquid water
namespace Water {

//! Specific heat capacity at constant pressure of water [J kg⁻¹ K⁻¹]
const double cp = 4186.84;

//! Latent heat of fusion of water/ice [J kg⁻¹]
const double Lf = Ice::Lf;

//! Latent heat of vaporization at 0˚C [J kg⁻¹]
const double Lv0 = Vapour::Lv0;

//! Proportionality constant between salinity in psu and freezing point depression [K psu⁻¹]
const double mu = 0.055;

//! Density of fresh water at 4˚C. [kg m⁻³]
const double rho = 1000;

//! Typical density of ocean water. [kg m⁻³]
const double rhoOcean = 1025;

//! Freezing point of pure water. [K]
const double Tf = Ice::Tm;

//! Freezing point of typical ocean water. [˚C]
const double TfOcean = -1.8;
}

namespace Nextsim {
//! Convert a temperature from ˚C to K
inline double kelvin(double celsius) { return celsius + PhysicalConstants::Tt; }

//! Convert a temperature from K to ˚C
inline double celsius(double kelvin) { return kelvin - PhysicalConstants::Tt; }

//! Convert an angle from radians to degrees
inline double degrees(double radians) { return radians * 360 / PhysicalConstants::tau; }

//! Convert an angle from degrees to radians
inline double radians(double degrees) { return degrees * PhysicalConstants::tau / 360; }

//! Convert a pressure from Pa to mbar
inline double mbar(double pascals) { return pascals / 100; }

//! Convert a pressure from mbar to Pa
inline double pascals(double mbar) { return mbar * 100; }
}
#endif /* SRC_INCLUDE_CONSTANTS_HPP */
