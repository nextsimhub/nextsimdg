/*!
 * @file PhysicsData.hpp
 * @date Sep 9, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef SRC_INCLUDE_PHYSICSDATA_HPP
#define SRC_INCLUDE_PHYSICSDATA_HPP

#include "BaseElementData.hpp"
#include "PrognosticData.hpp"

namespace Nextsim {

class PhysicsData : public BaseElementData {
public:
    PhysicsData() = default;
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
    //! Specific heat capacity of wet air
    inline double& heatCapacityWetAir() { return m_cspec; }
    //! Pressure due to wind drag [Pa]
    inline double& dragPressure() { return m_tau; }

    //! True ice thickness as updated [m]
    inline double& updatedIceTrueThickness() { return m_hi_new; }

    //! Mean thickness of snow (averaged over ice covered fraction) [m]
    inline double& updatedSnowTrueThickness() { return m_hs; }

    //! Updated value of the ice surface temperature [˚C]
    inline double& updatedIceSurfaceTemperature() { return m_TiceNew[0]; }

    inline double& updatedIceConcentration() { return m_conc_new; }

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
    std::array<double, N_ICE_TEMPERATURES> m_TiceNew;
    double m_conc_new; // updated ice concentration
};

} /* namespace Nextsim */

#endif /* SRC_INCLUDE_PHYSICSDATA_HPP */
