/*!
 * @file IThermodynamics.hpp
 *
 * @date Sep 29, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef SRC_INCLUDE_ITHERMODYNAMICS_HPP
#define SRC_INCLUDE_ITHERMODYNAMICS_HPP

namespace Nextsim {
class PrognosticData;
class PhysicsData;
class ExternalData;
class NextsimPhysics;

//! The interface class for ice thermodynamics.
class IThermodynamics {
public:
    virtual ~IThermodynamics() = default;

    /*!
     * @brief Calculate the ice thermodynamics.
     *
     * @param prog PrognosticData for this element (constant)
     * @param exter ExternalData for this element (constant)
     * @param phys PhysicsData for this element
     * @param nsphys Nextsim physics implementation data for this element.
     */
    virtual void calculate(const PrognosticData& prog, const ExternalData& exter, PhysicsData& phys,
        NextsimPhysics& nsphys)
        = 0;
};

}

#endif /* SRC_INCLUDE_ITHERMODYNAMICS_HPP */
