/*!
 * @file NextsimPhysics.hpp
 * @date Sep 9, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "BaseElementData.hpp"

#ifndef SRC_INCLUDE_NEXTSIMPHYSICS_HPP
#define SRC_INCLUDE_NEXTSIMPHYSICS_HPP

namespace Nextsim {

class PrognosticData;
class PhysicsData;
class ExternalData;

template<class Phys>
class ElementData;

class NextsimPhysics;

class NextsimPhysics : public BaseElementData {
public:
    NextsimPhysics() = default;
    ~NextsimPhysics() = default;

    inline void updateDerivedData(ElementData<NextsimPhysics>& data)
    {
        updateDerivedDataStatic(data);
    }
    inline void massFluxOpenWater(ElementData<NextsimPhysics>& data)
    {
        massFluxOpenWaterStatic(data);
    };
    inline void momentumFluxOpenWater(ElementData<NextsimPhysics>& data)
    {
        momentumFluxOpenWaterStatic(data);
    };
    inline void heatFluxOpenWater(ElementData<NextsimPhysics>& data)
    {
        heatFluxOpenWaterStatic(data);
    };
    inline void massFluxIceAtmosphere(ElementData<NextsimPhysics>& data)
    {
        massFluxIceAtmosphereStatic(data);
    };
    inline void heatFluxIceAtmosphere(ElementData<NextsimPhysics>& data)
    {
        heatFluxIceAtmosphereStatic(data);
    };
    inline void massFluxIceOcean(ElementData<NextsimPhysics>& data)
    {
        massFluxIceOceanStatic(data);
    };
    inline void heatFluxIceOcean(ElementData<NextsimPhysics>& data)
    {
        heatFluxIceOceanStatic(data);
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
    static void updateDerivedDataStatic(ElementData<NextsimPhysics>&);


    static void massFluxOpenWaterStatic(ElementData<NextsimPhysics>&);
    static void momentumFluxOpenWaterStatic(ElementData<NextsimPhysics>&);
    static void heatFluxOpenWaterStatic(ElementData<NextsimPhysics>&);

    static void massFluxIceAtmosphereStatic(ElementData<NextsimPhysics>&);
    static void heatFluxIceAtmosphereStatic(ElementData<NextsimPhysics>&);
    static void massFluxIceOceanStatic(ElementData<NextsimPhysics>&);
    static void heatFluxIceOceanStatic(ElementData<NextsimPhysics>&);
public:
    static double dragOcean_q;
    static double dragOcean_m(double windSpeed);
    static double dragOcean_t;

    static double latentHeatWater(double temperature);

    static SpecificHumidity specificHumidityWater;
    static SpecificHumidityIce specificHumidityIce;
};

} /* namespace Nextsim */

#endif /* SRC_INCLUDE_NEXTSIMPHYSICS_HPP */
