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
    static SpecificHumidity water();
    static SpecificHumidity ice();

    /*!
     * @brief Calculates humidity over fresh water or ice
     *
     * @param temperature Temperature of the water vapour [˚C]
     * @param pressure Hydrostatic pressure [Pa]
     */
    double operator()(const double temperature, const double pressure) const;
    /*!
     * @brief Calculates humidity over sea water.
     *
     * @param temperature Temperature of the water vapour [˚C]
     * @param pressure Hydrostatic pressure [Pa]
     * @param salinity Salinity of the liquid water [PSU]
     */
    double operator()(
        const double temperature, const double pressure, const double salinity) const;

    /*!
     * @brief Calculates humidity and its temperature dependence over fresh
     * water or ice.
     *
     * @param temperature Temperature of the water vapour [˚C]
     * @param pressure Hydrostatic pressure [Pa]
     */
    std::pair<double, double> valueAndDerivative(const double temperature, const double pressure) const;
    /*!
     * @brief Calculates humidity and its temperature dependence over sea
     * water.
     *
     * @param temperature Temperature of the water vapour [˚C]
     * @param pressure Hydrostatic pressure [Pa]
     * @param salinity Salinity of the liquid water [PSU]
     */
    std::pair<double, double> valueAndDerivative(
        const double temperature, const double pressure, const double salinity) const;

private:
    SpecificHumidity();
    const double m_a;
    const double m_b;
    const double m_c;
    const double m_d;
    const double m_bigA;
    const double m_bigB;
    const double m_bigC;
    const double m_alpha;
    const double m_beta;

};

} /* namespace Nextsim */

#endif /* SPECIFICHUMIDITY_HPP */
