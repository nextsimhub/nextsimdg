/*!
 * @file IPhysics1d.hpp
 *
 * @date Nov 23, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef SRC_INCLUDE_IPHYSICS1D_HPP
#define SRC_INCLUDE_IPHYSICS1D_HPP

namespace Nextsim {

class PrognosticData;
class ExternalData;
class PhysicsData;

//! Interface class for the column ice physics
class IPhysics1d {
    IPhysics1d() = 0;
    virtual ~IPhysics1d() = default;

    //! Perform the 1d physics calculation for this element, writing the data to the
    // PhysicsData argument.
    virtual void calculate(const PhysicsData&, const ExternalData&, PhysicsData&) = 0;
};
}
#endif /* SRC_INCLUDE_IPHYSICS1D_HPP */
