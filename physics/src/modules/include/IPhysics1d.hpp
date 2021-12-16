/*!
 * @file IPhysics1d.hpp
 *
 * @date Nov 23, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef SRC_INCLUDE_IPHYSICS1D_HPP
#define SRC_INCLUDE_IPHYSICS1D_HPP

#include "include/PhysicsData.hpp"
#include "include/PrognosticData.hpp"

namespace Nextsim {

class ExternalData;

//! The interface class for the column ice physics.
class IPhysics1d {
public:
    virtual ~IPhysics1d() = default;

    /*!
     * @brief Updates any derived quantities in PhysicsData.
     *
     * @details This function is declared virtual to be overridden if the implementing class needs to update
     * any class specific derived data.
     *
     * @param prog PrognosticData for this element (constant).
     * @param exter ExternalData for this element (constant).
     * @param phys PhysicsData for this element.
     */
    virtual void updateDerivedData(
        const PrognosticData& prog, const ExternalData& exter, PhysicsData& phys)
    {
        updateSpecificHumidityAir(exter, phys);
        updateSpecificHumidityWater(prog, exter, phys);
        updateSpecificHumidityIce(prog, exter, phys);

        updateAirDensity(exter, phys);
        updateHeatCapacityWetAir(exter, phys);

        phys.updatedSnowTrueThickness() = prog.snowTrueThickness();
        phys.updatedIceTrueThickness() = prog.iceTrueThickness();
    };

    /*!
     * @brief Performs the 1d physics calculation.
     *
     * @details Performs the one-dimensional physics calculation for this
     * element, writing the data to the PhysicsData argument.
     *
     * @param prog PrognosticData for this element (constant).
     * @param exter ExternalData for this element (constant).
     * @param phys PhysicsData for this element.
     */
    virtual void calculate(const PrognosticData&, const ExternalData&, PhysicsData&) = 0;

protected:
    /*!
     * @brief A virtual function that calculates the specific humidity in the
     * air.
     *
     * @param exter ExternalData for this element (constant).
     * @param phys PhysicsData for this element.
     */
    virtual void updateSpecificHumidityAir(const ExternalData& exter, PhysicsData& phys) = 0;
    /*!
     * @brief A virtual function that calculates the specific humidity at the
     * temperature of the sea surface and saturation.
     *
     * @param prog PrognosticData for this element (constant).
     * @param exter ExternalData for this element (constant).
     * @param phys PhysicsData for this element.
     */
    virtual void updateSpecificHumidityWater(
        const PrognosticData& prog, const ExternalData& exter, PhysicsData& phys)
        = 0;
    /*!
     * @brief A virtual function that calculates the specific humidity at the
     * temperature of the ice surface and saturation over ice.
     *
     * @param prog PrognosticData for this element (constant).
     * @param exter ExternalData for this element (constant).
     * @param phys PhysicsData for this element.
     */
    virtual void updateSpecificHumidityIce(
        const PrognosticData& prog, const ExternalData& exter, PhysicsData& phys)
        = 0;
    /*!
     * @brief A virtual function that calculates the air density.
     *
     * @param exter ExternalData for this element (constant).
     * @param phys PhysicsData for this element.
     */
    virtual void updateAirDensity(const ExternalData& exter, PhysicsData& phys) = 0;
    /*!
     * @brief A virtual function that calculates the heat capacity of the humid air.
     *
     * @param exter ExternalData for this element (constant).
     * @param phys PhysicsData for this element.
     */
    virtual void updateHeatCapacityWetAir(const ExternalData& exter, PhysicsData& phys) = 0;
};
}
#endif /* SRC_INCLUDE_IPHYSICS1D_HPP */
