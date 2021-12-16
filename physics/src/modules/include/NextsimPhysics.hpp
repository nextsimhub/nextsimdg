/*!
 * @file NextsimPhysics.hpp
 * @date Sep 9, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef SRC_INCLUDE_NEXTSIMPHYSICS_HPP
#define SRC_INCLUDE_NEXTSIMPHYSICS_HPP
#include <memory>

#include "include/BaseElementData.hpp"
#include "include/Configured.hpp"
#include "include/IPhysics1d.hpp"
#include "include/constants.hpp"

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
        OCEANALBEDO_KEY,
        I0_KEY,
        MINC_KEY,
        MINH_KEY,
    };

    void calculate(const PrognosticData&, const ExternalData&, PhysicsData&) override;

    //! Calculate the new ice formed this timestep on open water
    void newIceFormation(const PrognosticData& prog, const ExternalData& exter, PhysicsData& phys);
    //! The thickness of newly created ice in the current timestep
    inline double newIce() const { return m_newice; };

    //! Total ice-atmosphere heat flux [W m⁻²]
    inline double QIceAtmosphere() const { return m_Qia; };
    //! Ice to ocean heat flux [W m⁻²]
    inline double& QIceOceanHeat() { return m_Qio; };

    //! Rate of sublimation ice to vapour [kg s⁻¹ m⁻²]
    inline double sublimationRate() const { return m_subl; };

    //! Derivative of ice-atmosphere heat flux with respect to ice surface temperature [W m⁻² K⁻¹]
    inline double QDerivativeWRTTemperature() const { return m_dQ_dT; };

    //! Total amount of ice that resulted from the flooding of snow [m]
    inline double& totalIceFromSnow() { return m_hifroms; }

    //! Minimum ice concentration [1]
    static double minimumIceConcentration() { return minc; };
    //! Minimum ice thickness [m]
    static double minimumIceThickness() { return minh; };
    //! I0 parameter
    static double i0() { return m_I0; };

    // A class encapsulating the calculation of specific humidity.
    class SpecificHumidity {
    public:
        SpecificHumidity();
        /*!
         * @brief Calculates humidity over fresh water.
         *
         * @param temperature Temperature of the water vapour [˚C]
         * @param pressure Hydrostatic pressure [Pa]
         */
        double operator()(const double temperature, const double pressure) const;
        /*!
         * @brief Calculates humidity over sea water.
         *
         * @param temperature Temperature of the water vapour [˚C]
         * @param pressure Hydrostatic pressure [Pa]
         * @param salinity Salinity of the liquid water [PSU]
         */
        double operator()(
            const double temperature, const double pressure, const double salinity) const;

    protected:
        /*!
         * @brief Constructs an instance of the specific humidity calculator
         * with the specific parameters.
         *
         * @param a The 'a' scaling coefficient of the est factor.
         * @param b The 'b' exponential offset of the est factor.
         * @param c The 'c' temperature offset of the est factor.
         * @param d The 'd' temperature scaling of the est factor.
         * @param A The 'A' offset of the f factor.
         * @param B The 'B' pressure offset of the f factor.
         * @param C The 'C' quadratic temperature factor of the f factor.
         */
        SpecificHumidity(double a, double b, double c, double d, double A, double B, double C);
        /*!
         * @brief Calculates the f factor.
         *
         * @param temperature Water vapour temperature [˚C]
         * @param pressurePa Hydrostatic pressure [Pa]
         */
        double f(const double temperature, const double pressurePa) const;
        /*!
         * @brief Calculates the est factor.
         *
         * @param temperature Water vapour temperature [˚C]
         * @param salinity Liquid water salinity [PSU]
         */
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
    // A class encapsulating the calculation of specific humidity over ice.
    class SpecificHumidityIce : public SpecificHumidity {
    public:
        SpecificHumidityIce();
        /*!
         * @brief Calculates humidity over ice.
         *
         * @param temperature Temperature of the water vapour [˚C]
         * @param pressure Hydrostatic pressure [Pa]
         */
        double operator()(const double temperature, const double pressure) const;
        /*!
         * @brief Derivative of the specific humdity over ice with respect to
         * temperature.
         *
         * @param temperature Temperature of the water vapour [˚C]
         * @param pressure Hydrostatic pressure [Pa]
         */
        double dq_dT(const double temperature, const double pressure) const;
    };

protected:
    /*!
     * @brief Calculates the specific humidity in the air.
     *
     * @param exter ExternalData for this element (constant).
     * @param phys PhysicsData for this element.
     */
    void updateSpecificHumidityAir(const ExternalData& exter, PhysicsData& phys) override;
    /*!
     * @brief Calculates the specific humidity at the temperature of the sea
     * surface and saturation.
     *
     * @param prog PrognosticData for this element (constant).
     * @param exter ExternalData for this element (constant).
     * @param phys PhysicsData for this element.
     */
    void updateSpecificHumidityWater(
        const PrognosticData& prog, const ExternalData& exter, PhysicsData& phys) override;
    /*!
     * @brief Calculates the specific humidity at the temperature of the ice
     * surface and saturation over ice.
     *
     * @param prog PrognosticData for this element (constant).
     * @param exter ExternalData for this element (constant).
     * @param phys PhysicsData for this element.
     */
    void updateSpecificHumidityIce(
        const PrognosticData& prog, const ExternalData& exter, PhysicsData& phys) override;
    /*!
     * @brief Calculates the air density.
     *
     * @param exter ExternalData for this element (constant).
     * @param phys PhysicsData for this element.
     */
    void updateAirDensity(const ExternalData& exter, PhysicsData& phys) override;
    /*!
     * @brief Calculates the heat capacity of the humid air.
     *
     * @param exter ExternalData for this element (constant).
     * @param phys PhysicsData for this element.
     */
    void updateHeatCapacityWetAir(const ExternalData& exter, PhysicsData& phys) override;

private:
    void massFluxOpenWater(PhysicsData& phys);
    void momentumFluxOpenWater(PhysicsData& phys);
    void heatFluxOpenWater(
        const PrognosticData& prog, const ExternalData& exter, PhysicsData& phys);

    void massFluxIceAtmosphere(const PrognosticData& prog, PhysicsData& phys);
    void heatFluxIceAtmosphere(
        const PrognosticData& prog, const ExternalData& exter, PhysicsData& phys);
    void massFluxIceOcean(const PrognosticData& prog, const ExternalData& exter, PhysicsData& phys);
    void heatFluxIceOcean(const PrognosticData& prog, const ExternalData& exter, PhysicsData& phys);
    void lateralGrowth(const PrognosticData& prog, const ExternalData& exter, PhysicsData& phys);

    // Phase change rates
    double m_evap;
    double m_subl;

    // Open water heat fluxes
    double m_Qow;
    double m_Qlwow;
    double m_Qswow;
    double m_Qlhow;
    double m_Qshow;

    // Ice heat fluxes
    double m_Qi;
    double m_Qlwi;
    double m_Qswi;
    double m_Qlhi;
    double m_Qshi;
    double m_dQ_dT;

    // ice-ocean fluxes
    double m_Qio;

    //! ice-atmosphere fluxes
    double m_Qia;

    // Thickness of ice generated from flooding of snow [m]
    double m_hifroms;

    static double dragOcean_q;
    static double dragOcean_m(double windSpeed);
    static double dragOcean_t;
    static double dragIce_t;

    static double m_oceanAlbedo;

    // Ice-ocean heat flux
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
