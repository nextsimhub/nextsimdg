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
class UnusedData;

template<class Phys>
class ElementData;

class NextsimPhysics;

class NextsimPhysics : public BaseElementData {
public:
    NextsimPhysics() = default;
    ~NextsimPhysics() = default;

    inline static void updateDerivedData(
            const PrognosticData& prog,
            const ExternalData& exter,
            PhysicsData& phys,
            const UnusedData&)
    {
        updateDerivedDataStatic(prog, exter, phys);
    }
    inline static void massFluxOpenWater(
            const UnusedData&,
            const UnusedData&,
            PhysicsData& phys,
            const UnusedData&)
    {
        massFluxOpenWaterStatic(phys);
    };
    inline static void momentumFluxOpenWater(
            const UnusedData&,
            const UnusedData&,
            PhysicsData& phys,
            const UnusedData&)
    {
        momentumFluxOpenWaterStatic(phys);
    };
    inline static void heatFluxOpenWater(
            const PrognosticData& prog,
            const ExternalData& exter,
            PhysicsData& phys,
            UnusedData&)
    {
        heatFluxOpenWaterStatic(prog, exter, phys);
    }
    inline static void massFluxIceAtmosphere(
            const PrognosticData& prog,
            const ExternalData& exter,
            PhysicsData& phys,
            UnusedData&)
    {
        massFluxIceAtmosphereStatic(prog, exter, phys);
    };
    inline static void heatFluxIceAtmosphere(
            const PrognosticData& prog,
            const ExternalData& exter,
            PhysicsData& phys,
            UnusedData&)
    {
        heatFluxIceAtmosphereStatic(prog, exter, phys);
    };
    inline static void massFluxIceOcean(
            const PrognosticData& prog,
            const ExternalData& exter,
            PhysicsData& phys,
            UnusedData&)
    {
        massFluxIceOceanStatic(prog, exter, phys);
    };
    inline static void heatFluxIceOcean(
            const PrognosticData& prog,
            const ExternalData& exter,
            PhysicsData& phys,
            UnusedData&)
    {
        heatFluxIceOceanStatic(prog, exter, phys);
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
    static void updateDerivedDataStatic(
            const PrognosticData& prog,
            const ExternalData& exter,
            PhysicsData& phys);

    static void massFluxOpenWaterStatic(PhysicsData& phys);
    static void momentumFluxOpenWaterStatic(PhysicsData& phys);
    static void heatFluxOpenWaterStatic(
            const PrognosticData& prog,
            const ExternalData& exter,
            PhysicsData& phys);

    static void massFluxIceAtmosphereStatic(
            const PrognosticData& prog,
            const ExternalData& exter,
            PhysicsData& phys);
    static void heatFluxIceAtmosphereStatic(
            const PrognosticData& prog,
            const ExternalData& exter,
            PhysicsData& phys);
    static void massFluxIceOceanStatic(
            const PrognosticData& prog,
            const ExternalData& exter,
            PhysicsData& phys);
    static void heatFluxIceOceanStatic(
            const PrognosticData& prog,
            const ExternalData& exter,
            PhysicsData& phys);
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
