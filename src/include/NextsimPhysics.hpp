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

class ElementData;

class NextsimPhysics : public BaseElementData {
public:
    NextsimPhysics() = default;
    ~NextsimPhysics() = default;

    inline void updateDerivedData(ElementData& data)
    {
        updateDerivedDataStatic(data);
    }
    inline void massFluxOpenWater(ElementData& data)
    {
        massFluxOpenWaterStatic(data);//
    };
    inline void momentumFluxOpenWater(ElementData& data)
    {
        momentumFluxOpenWaterStatic(data);//(data, data);
    };
    inline void heatFluxOpenWater(ElementData& data)
    {
        heatFluxOpenWaterStatic(data);//(data, data, data);
    };
    inline void massFluxIceAtmosphere(ElementData& data)
    {
        massFluxIceAtmosphereStatic(data);//(data, data, data);
    };
    inline void heatFluxIceAtmosphere(ElementData& data)
    {
        heatFluxIceAtmosphereStatic(data);//(data, data, data);
    };
    inline void massFluxIceOcean(ElementData& data)
    {
        massFluxIceOceanStatic(data);//(data, data, data);
    };
    inline void heatFluxIceOcean(ElementData& data)
    {
        heatFluxIceOceanStatic(data);//(data, data, data);
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
    static void updateDerivedDataStatic(ElementData&);
    //static void updateDerivedData(const PrognosticData& prog, const ExternalData &exter, PhysicsData& phys);

    static void massFluxOpenWaterStatic(ElementData&);//const PrognosticData& prog, PhysicsData& phys);
    static void momentumFluxOpenWaterStatic(ElementData&);//(const PrognosticData& prog, PhysicsData& phys);
    static void heatFluxOpenWaterStatic(ElementData&);//(const PrognosticData& prog, const ExternalData &exter, PhysicsData& phys);

    static void massFluxIceAtmosphereStatic(ElementData&);//(const PrognosticData& prog, const ExternalData &exter, PhysicsData& phys);
    static void heatFluxIceAtmosphereStatic(ElementData&);//(const PrognosticData& prog, const ExternalData &exter, PhysicsData& phys);
    static void massFluxIceOceanStatic(ElementData&);//(const PrognosticData& prog, const ExternalData &exter, PhysicsData& phys);
    static void heatFluxIceOceanStatic(ElementData&);//(const PrognosticData& prog, const ExternalData &exter, PhysicsData& phys);
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
