/*!
 * @file ElementData.hpp
 * @date Sep 9, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef SRC_INCLUDE_ELEMENTDATA_HPP
#define SRC_INCLUDE_ELEMENTDATA_HPP

#include "PhysicsData.hpp"
#include "PrognosticData.hpp"
#include "ExternalData.hpp"

#if 1 // TODO: A more sophisticated system of selecting class that implements the physics.
#include "NextsimPhysics.hpp"
typedef Nextsim::NextsimPhysics PhysicsImpl;
#else
namespace Nextsim {
class NoPhysics: public BaseElementData {

};
typedef NoPhysics PhysicsImpl;
}
#endif

namespace Nextsim {

class ElementData: public PrognosticData,
                   public PhysicsData,
                   public ExternalData,
                   public PhysicsImpl {
public:
    ElementData() = default;
    ~ElementData() = default;

    using PrognosticData::operator=;
    using PhysicsData::operator=;
    using ExternalData::operator=;
};

} /* namespace Nextsim */

#endif /* SRC_INCLUDE_ELEMENTDATA_HPP */
