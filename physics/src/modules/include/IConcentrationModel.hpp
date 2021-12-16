/*!
 * @file IConcentrationModel.hpp
 *
 * @date Nov 11, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef SRC_INCLUDE_ICONCENTRATIONMODEL_HPP
#define SRC_INCLUDE_ICONCENTRATIONMODEL_HPP

namespace Nextsim {

class PhysicsData;
class PrognosticData;
class NextsimPhysics;

//! The interface class for ice concentration update calculations.
class IConcentrationModel {
public:
    virtual ~IConcentrationModel() = default;

    /*!
     * @brief A virtual function that calculates the amount of freezing during
     * the timestep.
     *
     * @param prog PrognosticData for this element (constant).
     * @param phys PhysicsData for this element.
     * @param nsphys Nextsim physics implementation data for this element.
     */
    virtual double freeze(const PrognosticData&, PhysicsData&, NextsimPhysics&) const = 0;
    /*!
     * @brief A virtual function that calculates the amount of melting during
     * the timestep.
     *
     * @param prog PrognosticData for this element (constant).
     * @param phys PhysicsData for this element.
     * @param nsphys Nextsim physics implementation data for this element.
     */
    virtual double melt(const PrognosticData&, PhysicsData&, NextsimPhysics&) const = 0;
};
}

#endif /* SRC_INCLUDE_ICONCENTRATIONMODEL_HPP */
