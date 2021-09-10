/*!
 * @file NextsimPhysics.hpp
 * @date Sep 9, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef SRC_INCLUDE_NEXTSIMPHYSICS_HPP
#define SRC_INCLUDE_NEXTSIMPHYSICS_HPP

namespace Nextsim {

class ElementData;
class PrognosticData;
class PhysicsData;

class NextsimPhysics {
public:
    NextsimPhysics();
    virtual ~NextsimPhysics();

    inline void updateDerivedData(ElementData& data)
    {
        updateDerivedData(data, data);
    }
    inline void massFluxOpenWater(ElementData& data)
    {
        massFluxOpenWater(data, data);
    };
    inline void momentumFluxOpenWater(ElementData& data)
    {
        momentumFluxOpenWater(data, data);
    };
    inline void heatFluxOpenWater(ElementData& data)
    {
        heatFluxOpenWater(data, data);
    };

private:
    static void updateDerivedData(const PrognosticData& prog, PhysicsData& phys);
    static void massFluxOpenWater(const PrognosticData& prog, PhysicsData& phys);
    static void momentumFluxOpenWater(const PrognosticData& prog, PhysicsData& phys);
    static void heatFluxOpenWater(const PrognosticData& prog, PhysicsData& phys);
};

} /* namespace Nextsim */

#endif /* SRC_INCLUDE_NEXTSIMPHYSICS_HPP */
