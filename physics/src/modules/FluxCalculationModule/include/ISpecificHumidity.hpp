/*
 *
 * @author  Tim Spain <timothy.spain@nersc.no>
 */

#ifndef ISPECIFICHUMIDITY_HPP
#define ISPECIFICHUMIDITY_HPP

#include "include/FloatType.hpp"

#include <utility>

namespace Nextsim {
//! An interface class for the specific humidity calculations.
class ISpecificHumidity {
public:
    virtual ~ISpecificHumidity() = default;

    /*!
     * @brief Calculates humidity over fresh water or ice
     *
     * @param temperature Temperature of the water vapour [˚C]
     * @param pressure Hydrostatic pressure [Pa]
     */
    virtual FloatType operator()(FloatType temperature, FloatType pressure) const = 0;
    /*!
     * @brief Calculates humidity over sea water.
     *
     * @param temperature Temperature of the water vapour [˚C]
     * @param pressure Hydrostatic pressure [Pa]
     * @param salinity Salinity of the liquid water [PSU]
     */
    virtual FloatType operator()(
        FloatType temperature, FloatType pressure, FloatType salinity) const
        = 0;

    /*!
     * @brief Calculates humidity and its temperature dependence over fresh
     * water or ice.
     *
     * @param temperature Temperature of the water vapour [˚C]
     * @param pressure Hydrostatic pressure [Pa]
     */
    virtual std::pair<FloatType, FloatType> valueAndDerivative(
        FloatType temperature, FloatType pressure) const
        = 0;
    /*!
     * @brief Calculates humidity and its temperature dependence over sea
     * water.
     *
     * @param temperature Temperature of the water vapour [˚C]
     * @param pressure Hydrostatic pressure [Pa]
     * @param salinity Salinity of the liquid water [PSU]
     */
    virtual std::pair<FloatType, FloatType> valueAndDerivative(
        FloatType temperature, FloatType pressure, FloatType salinity) const
        = 0;
};

} /* namespace Nextsim */

#endif /* ISPECIFICHUMIDITY_HPP */
