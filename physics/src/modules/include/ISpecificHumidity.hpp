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
    // Default constructor is for water
    ISpecificHumidity()
    : waterPtr(nullptr)
    , icePtr(nullptr)
    {
    }
    virtual ~ISpecificHumidity() = default;

    //! @brief Returns a reference to the static water instance
    static ISpecificHumidity& water() { return *waterPtr; }
    //! @brief Returns a reference to the static ice instance
    static ISpecificHumidity& ice() { return *icePtr; }

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
    virtual std::pair<double, double> valueAndDerivative(double temperature, double pressure) const = 0;
    /*!
     * @brief Calculates humidity and its temperature dependence over sea
     * water.
     *
     * @param temperature Temperature of the water vapour [˚C]
     * @param pressure Hydrostatic pressure [Pa]
     * @param salinity Salinity of the liquid water [PSU]
     */
    virtual std::pair<double, double> valueAndDerivative(double temperature,
            double pressure, double salinity) const = 0;

protected:
    ISpecificHumidity* waterPtr;
    ISpecificHumidity* icePtr;
};

} /* namespace Nextsim */

#endif /* ISPECIFICHUMIDITY_HPP */
