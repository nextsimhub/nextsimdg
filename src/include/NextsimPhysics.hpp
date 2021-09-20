/*!
 * @file NextsimPhysics.hpp
 * @date Sep 9, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef SRC_INCLUDE_NEXTSIMPHYSICS_HPP
#define SRC_INCLUDE_NEXTSIMPHYSICS_HPP

#include "ElementData.hpp"

namespace Nextsim {

class PrognosticData;
class ExternalData;
class PhysicsData;

class NextsimPhysics {
public:
    NextsimPhysics();
    virtual ~NextsimPhysics();

    inline void updateDerivedData(ElementData& data)
    {
        updateDerivedData(data, data, data);
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
    static void setDragOcean_q(double dragOcean_q);
    static void setDragOcean_t(double dragOcean_t);

    class SpecificHumidity {
    public:
        SpecificHumidity();
        double operator()(const double temperature, const double pressure) const;
        double operator()(const double temperature, const double pressure, const double salinity) const;
    protected:
        SpecificHumidity(double, double, double, double, double, double, double);
        double f(const double pressurePa, const double temperature) const;
        double est(const double temperature, const double salinity) const;
        const double m_a;
        const double m_b;
        const double m_c;
        const double m_d;
        const double m_bigA;
        const double m_bigB;
        const double m_bigC;
        const double m_alpha;
        const double m_beta;
    };
    class SpecificHumidityIce : public SpecificHumidity {
    public:
           SpecificHumidityIce();
    protected:
           double operator()(const double temperature, const double pressure) const;
           double dq_dT(const double temperature, const double pressure) const;
    };
private:
    static void updateDerivedData(const PrognosticData& prog, const ExternalData &exter, PhysicsData& phys);
    static void massFluxOpenWater(const PrognosticData& prog, PhysicsData& phys);
    static void momentumFluxOpenWater(const PrognosticData& prog, PhysicsData& phys);
    static void heatFluxOpenWater(const PrognosticData& prog, const ExternalData &exter, PhysicsData& phys);

    static double dragOcean_q;
    static double dragOcean_m(double windSpeed);
    static double dragOcean_t;

    static double latentHeatWater(double temperature);
    static SpecificHumidity specificHumidityWater;

    static SpecificHumidityIce specificHumidityIce;
};

} /* namespace Nextsim */

#endif /* SRC_INCLUDE_NEXTSIMPHYSICS_HPP */
