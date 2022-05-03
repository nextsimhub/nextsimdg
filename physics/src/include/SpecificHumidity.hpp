/*
 * @file SpecificHumidity.hpp
 *
 * @date May 2, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef SPECIFICHUMIDITY_HPP
#define SPECIFICHUMIDITY_HPP

#include <utility>

namespace Nextsim {

class SpecificHumidity {
public:
    //! @brief Returns a reference to the static water instance
    static SpecificHumidity& water();
    //! @brief Returns a reference to the static ice instance
    static SpecificHumidity& ice();

    /*!
     * @brief Calculates humidity over fresh water or ice
     *
     * @param temperature Temperature of the water vapour [˚C]
     * @param pressure Hydrostatic pressure [Pa]
     */
    double operator()(double temperature, double pressure) const;
    /*!
     * @brief Calculates humidity over sea water.
     *
     * @param temperature Temperature of the water vapour [˚C]
     * @param pressure Hydrostatic pressure [Pa]
     * @param salinity Salinity of the liquid water [PSU]
     */
    double operator()(double temperature, double pressure, double salinity) const;

    /*!
     * @brief Calculates humidity and its temperature dependence over fresh
     * water or ice.
     *
     * @param temperature Temperature of the water vapour [˚C]
     * @param pressure Hydrostatic pressure [Pa]
     */
    std::pair<double, double> valueAndDerivative(double temperature, double pressure) const;
    /*!
     * @brief Calculates humidity and its temperature dependence over sea
     * water.
     *
     * @param temperature Temperature of the water vapour [˚C]
     * @param pressure Hydrostatic pressure [Pa]
     * @param salinity Salinity of the liquid water [PSU]
     */
    std::pair<double, double> valueAndDerivative(double temperature,
            double pressure, double salinity) const;

private:
    // Default constructor is for water
    SpecificHumidity();
    // General constructor
    SpecificHumidity(double a, double b, double c, double d, double bigA, double bigB, double bigC);

    /*!
     * @brief Performs the calculation
     *
     * @param temperature Temperature of the water vapour [˚C]
     * @param pressure Hydrostatic pressure [Pa]
     * @param salinity Salinity of the liquid water [PSU]
     * @param doDeriv perform the derivative calculation
     */
    std::pair<double, double> calculate(double temperature, double pressure, double salinity, bool doDeriv) const;

    /*!
     * @brief Calculates the f factor.
     *
     * @param temperature Water vapour temperature [˚C]
     * @param pressurePa Hydrostatic pressure [Pa]
     */
    double f(double temperature, double pressurePa) const;
    /*!
     * @brief Calculates the est factor.
     *
     * @param temperature Water vapour temperature [˚C]
     * @param salinity Liquid water salinity [PSU]
     */
    double est(double temperature, double salinity) const;

    const double m_a;
    const double m_b;
    const double m_c;
    const double m_d;
    const double m_bigA;
    const double m_bigB;
    const double m_bigC;
    const double m_alpha;
    const double m_beta;

    static SpecificHumidity waterInst;
    static SpecificHumidity iceInst;
};

} /* namespace Nextsim */

#endif /* SPECIFICHUMIDITY_HPP */
