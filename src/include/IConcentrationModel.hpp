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

class IConcentrationModel {
public:
    virtual ~IConcentrationModel() = default;

    virtual double freeze(const PrognosticData&, NextsimPhysics&) const = 0;
    virtual double melt(const PrognosticData&, NextsimPhysics&) const = 0;
};
}

#endif /* SRC_INCLUDE_ICONCENTRATIONMODEL_HPP */
