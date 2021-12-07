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

//! Interface class for the column ice physics
class IPhysics1d {
public:
    virtual ~IPhysics1d() = default;

    /*!
     * Update any derived quantities in PhysicsData.
     *
     * This function is declared virtual to be overridden if the implementing class needs to update
     * any class specific derived data.
     *
     * @param prog PrognosticData for this cell (constant)
     * @param exter ExternalData for this cell (constant)
     * @param phys PhysicsData for this point
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

    //! Perform the 1d physics calculation for this element, writing the data to the
    // PhysicsData argument.
    virtual void calculate(const PrognosticData&, const ExternalData&, PhysicsData&) = 0;

protected:
    virtual void updateSpecificHumidityAir(const ExternalData& exter, PhysicsData& phys) = 0;
    virtual void updateSpecificHumidityWater(
        const PrognosticData& prog, const ExternalData& exter, PhysicsData& phys)
        = 0;
    virtual void updateSpecificHumidityIce(
        const PrognosticData& prog, const ExternalData& exter, PhysicsData& phys)
        = 0;
    virtual void updateAirDensity(const ExternalData& exter, PhysicsData& phys) = 0;
    virtual void updateHeatCapacityWetAir(const ExternalData& exter, PhysicsData& phys) = 0;
};
}
#endif /* SRC_INCLUDE_IPHYSICS1D_HPP */
