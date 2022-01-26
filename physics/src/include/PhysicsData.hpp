/*!
 * @file PhysicsData.hpp
 * @date Sep 9, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef SRC_INCLUDE_PHYSICSDATA_HPP
#define SRC_INCLUDE_PHYSICSDATA_HPP

#include "include/BaseElementData.hpp"
#include "include/IPrognosticUpdater.hpp"
#include "include/PrognosticData.hpp"

#include <vector>
namespace Nextsim {

//! A class holding common physics data.
class PhysicsData : public BaseElementData, public IPrognosticUpdater {
public:
    PhysicsData()
        : PhysicsData(1)
    {
    }
    PhysicsData(int nIceLayers)
        : m_rho(0)
        , m_wspeed(0)
        , m_sphumi(0)
        , m_sphuma(0)
        , m_sphumw(0)
        , m_cspec(0)
        , m_tau(0)
        , m_conc_new(0)
        , m_hi_new(0)
        , m_hs(0)
        , m_TiceNew(nIceLayers, 0.)
    {
    }

    ~PhysicsData() = default;

    //! Density of air at the current temperature and humidity [kg m⁻³]
    inline double& airDensity() { return m_rho; };
    //! Wind speed [m s⁻¹]
    inline double& windSpeed() { return m_wspeed; }
    //! Specific humidity over the water [kg kg⁻¹]
    inline double& specificHumidityWater() { return m_sphumw; }
    //! Specific humidity over the ice [kg kg⁻¹]
    inline double& specificHumidityIce() { return m_sphumi; }
    //! Specific humidity of the air [kg kg⁻¹]
    inline double& specificHumidityAir() { return m_sphuma; }
    //! Mixing ratio of water vapour in the air [kg kg⁻¹]
    inline double mixingRatio() { return m_sphuma / (1 - m_sphuma); }
    //! Specific heat capacity of wet air [J kg⁻¹ K⁻¹]
    inline double& heatCapacityWetAir() { return m_cspec; }
    //! Pressure due to wind drag [Pa]
    inline double& dragPressure() { return m_tau; }

    //! True ice thickness as updated [m]
    inline double& updatedIceTrueThickness() { return m_hi_new; }
    //! Mean ice thickness, as updated [m]
    double updatedIceThickness() const override { return m_hi_new * m_conc_new; }

    //! Mean thickness of snow (averaged over ice covered fraction) [m]
    inline double& updatedSnowTrueThickness() { return m_hs; }
    //! Mean thickness of snow (averaged over data element) [m]
    double updatedSnowThickness() const override { return m_hs * m_conc_new; }

    //! Updated value of the ice surface temperature [˚C]
    inline double& updatedIceSurfaceTemperature() { return m_TiceNew[0]; }
    //! Updated layer-wise ice temperatures [˚C]
    const std::vector<double>& updatedIceTemperatures() const override { return m_TiceNew; };

    //! Updated value of the ice concentration [1]
    inline double& updatedIceConcentration() { return m_conc_new; }
    //! Updated value of the ice concentration [1]
    double updatedIceConcentration() const override { return m_conc_new; }

private:
    double m_rho;
    double m_wspeed;
    double m_sphumw;
    double m_sphumi;
    double m_sphuma;
    double m_cspec;
    double m_tau;

    // thermodynamic values
    double m_hi_new; // updated true ice thickness [m]
    double m_hs;
    std::vector<double> m_TiceNew;
    double m_conc_new; // updated ice concentration
};

} /* namespace Nextsim */

#endif /* SRC_INCLUDE_PHYSICSDATA_HPP */
