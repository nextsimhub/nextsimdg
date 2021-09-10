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

    static class SpecificHumidity {
    public:
        SpecificHumidity();
        double operator()(double temperature, double pressure);
        double operator()(double temperature, double pressure, double salinity);
    protected:
        double m_a;
        double m_b;
        double m_c;
        double m_d;
        double m_bigA;
        double m_bigB;
        double m_bigC;
        double alpha;
        double beta;
    } specificHumidityWater;

    static class SpecificHumidityIce : public SpecificHumidity {
        SpecificHumidityIce();
        double operator()(double temperature, double pressure);
    } specificHumidityIce;
};

} /* namespace Nextsim */

#endif /* SRC_INCLUDE_NEXTSIMPHYSICS_HPP */
