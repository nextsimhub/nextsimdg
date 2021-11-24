/*!
 * @file NextsimPhysics.hpp
 * @date Sep 9, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef SRC_INCLUDE_NEXTSIMPHYSICS_HPP
#define SRC_INCLUDE_NEXTSIMPHYSICS_HPP
#include <memory>

#include "BaseElementData.hpp"
#include "Configured.hpp"
#include "constants.hpp"
#include "IPhysics1d.hpp"

namespace Nextsim {

class PrognosticData;
class PhysicsData;
class ExternalData;
class UnusedData;

class IIceAlbedo;
class IThermodynamics;
class IIceOceanHeatFlux;
class IConcentrationModel;

template <class Phys> class ElementData;

class NextsimPhysics;

class NextsimPhysics : public BaseElementData,
                       public Configured<NextsimPhysics>,
                       public IPhysics1d {
public:
    NextsimPhysics();

    void configure() override;
    enum {
        DRAGOCEANQ_KEY,
        DRAGOCEANT_KEY,
        DRAGICET_KEY,
        I0_KEY,
        MINC_KEY,
        MINH_KEY,
    };

    void calculate(const PrognosticData&, const ExternalData&, PhysicsData&) override;

    //! Calculate the new ice formed this timestep on open water
    void newIceFormation(const PrognosticData& prog, const ExternalData& exter, PhysicsData& phys);
    //! The thickness of newly created ice in the current timestep
    inline double newIce() const { return m_newice; };


    static double minimumIceConcentration() { return minc; };
    static double minimumIceThickness() { return minh; };
    static double i0() { return m_I0; };

    class SpecificHumidity {
    public:
        SpecificHumidity();
        double operator()(const double temperature, const double pressure) const;
        double operator()(
            const double temperature, const double pressure, const double salinity) const;

    protected:
        SpecificHumidity(double, double, double, double, double, double, double);
        double f(const double temperature, const double pressurePa) const;
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
        double operator()(const double temperature, const double pressure) const;
        double dq_dT(const double temperature, const double pressure) const;
    };
protected:
    void updateSpecificHumidityAir(const ExternalData& exter, PhysicsData& phys) override;
    void updateSpecificHumidityWater(
            const PrognosticData& prog, const ExternalData& exter, PhysicsData& phys) override;
    void updateSpecificHumidityIce(
            const PrognosticData& prog, const ExternalData& exter, PhysicsData& phys) override;
    void updateAirDensity(const ExternalData& exter, PhysicsData& phys) override;
    void updateHeatCapacityWetAir(const ExternalData& exter, PhysicsData& phys) override;

private:
    void massFluxOpenWater(PhysicsData& phys);
    void momentumFluxOpenWater(PhysicsData& phys);
    void heatFluxOpenWater(
        const PrognosticData& prog, const ExternalData& exter, PhysicsData& phys);

    void massFluxIceAtmosphere(const PrognosticData& prog, PhysicsData& phys);
    void heatFluxIceAtmosphere(
        const PrognosticData& prog, const ExternalData& exter, PhysicsData& phys);
    void massFluxIceOcean(const PrognosticData& prog, const ExternalData& exter,
        PhysicsData& phys);
    void heatFluxIceOcean(const PrognosticData& prog, const ExternalData& exter, PhysicsData& phys);
    void lateralGrowth(const PrognosticData& prog, const ExternalData& exter, PhysicsData& phys);

    static double dragOcean_q;
    static double dragOcean_m(double windSpeed);
    static double dragOcean_t;
    static double dragIce_t;

    // Private instance functions

    // Ice-ocean heat flux
    double m_Qio; // Ice-ocean heat flux
    static IIceOceanHeatFlux* iceOceanHeatFluxImpl;
    // New ice created by cooling below freezing
    double m_newice;

    static IConcentrationModel* iConcentrationModelImpl;

    static double m_I0;

    static double latentHeatWater(double temperature);
    static double latentHeatIce(double temperature);

    static double minc; // minimum ice concentration
    static double minh; // minimum ice true thickness [m]

    static SpecificHumidity specHumWater;
    static SpecificHumidityIce specHumIce;

    static IIceAlbedo* iIceAlbedoImpl;
    static IThermodynamics* iThermo;
};

} /* namespace Nextsim */

#endif /* SRC_INCLUDE_NEXTSIMPHYSICS_HPP */
