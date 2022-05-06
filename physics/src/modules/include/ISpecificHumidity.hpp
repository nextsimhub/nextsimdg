/*
 * @file ISpecificHumidity.hpp
 *
 * @date May 2, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef ISPECIFICHUMIDITY_HPP
#define ISPECIFICHUMIDITY_HPP

#include <utility>

namespace Nextsim {

class ISpecificHumidity {
public:
    virtual ~ISpecificHumidity() = default;

    /*!
     * @brief Calculates humidity over fresh water or ice
     *
     * @param temperature Temperature of the water vapour [˚C]
     * @param pressure Hydrostatic pressure [Pa]
     */
    virtual double operator()(double temperature, double pressure) const = 0;
    /*!
     * @brief Calculates humidity over sea water.
     *
     * @param temperature Temperature of the water vapour [˚C]
     * @param pressure Hydrostatic pressure [Pa]
     * @param salinity Salinity of the liquid water [PSU]
     */
    virtual double operator()(double temperature, double pressure, double salinity) const = 0;

    /*!
     * @brief Calculates humidity and its temperature dependence over fresh
     * water or ice.
     *
     * @param temperature Temperature of the water vapour [˚C]
     * @param pressure Hydrostatic pressure [Pa]
     */
    virtual std::pair<double, double> valueAndDerivative(
        double temperature, double pressure) const = 0;
    /*!
     * @brief Calculates humidity and its temperature dependence over sea
     * water.
     *
     * @param temperature Temperature of the water vapour [˚C]
     * @param pressure Hydrostatic pressure [Pa]
     * @param salinity Salinity of the liquid water [PSU]
     */
    virtual std::pair<double, double> valueAndDerivative(
        double temperature, double pressure, double salinity) const = 0;
};

} /* namespace Nextsim */

#endif /* ISPECIFICHUMIDITY_HPP */
