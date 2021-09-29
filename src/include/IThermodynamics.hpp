/*
 * @file IThermodynamics.hpp
 *
 * @date Sep 29, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef SRC_INCLUDE_ITHERMODYNAMICS_HPP
#define SRC_INCLUDE_ITHERMODYNAMICS_HPP

class PrognosticData;
class PhysicsData;
class ExternalData;
class NextsimPhysics;

class IThermodynamics {
public:
    virtual ~IThermodynamics() = default;

    virtual void calculate(
            const PrognosticData& prog,
            const ExternalData& exter,
            PhysicsData& phys,
            NextsimPhysics& nsphys);
};



#endif /* SRC_INCLUDE_ITHERMODYNAMICS_HPP */
