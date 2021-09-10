/*!
 * @file PhysicsData.hpp
 * @date Sep 9, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef SRC_INCLUDE_PHYSICSDATA_HPP
#define SRC_INCLUDE_PHYSICSDATA_HPP

#include <BaseElementData.hpp>

namespace Nextsim {

class PhysicsData: public BaseElementData {
public:
    PhysicsData();
    virtual ~PhysicsData();

    //! Rate of evaporation (+ve) or condensation (-ve) [kg s⁻¹ m⁻²]
    inline double& evaporationRate()
    {
        return m_evap;
    };
    //! Density of air at the current temperature and humidity [kg m⁻³]
    inline double& airDensity()
    {
        return m_rho;
    };
    //! Wind speed [m s⁻¹]
    inline double& windSpeed()
    {
        return m_wspeed;
    }
    //! Specific humidity over the water [kg kg⁻¹]
    inline double& specificHumidityWater()
    {
        return m_sphumw;
    }
    //! Specific humidity of the air [kg kg⁻¹]
    inline double& specificHumidityAir()
    {
        return m_sphuma;
    }
    //! Mixing ratio of water vapour in the air [kg kg⁻¹]
    inline double mixingRatio( )
    {
        return m_sphuma / (1 - m_sphuma);
    }
    //! Pressure due to wind drag [Pa]
    inline double& dragPressure()
    {
        return m_tau;
    }
    //! Ocean drag scaling factor for humidity [s] (apparently)
    static inline double& dragOcean_q()
    {
        return m_drag_ocean_q;
    }
    //! Specific drag from density to pressure [m s⁻²]
    double dragOcean_m();
private:
    double m_evap;
    double m_rho;
    double m_wspeed;
    double m_sphumw;
    double m_sphuma;
    double m_tau;

    static double m_drag_ocean_q;
};

} /* namespace Nextsim */

#endif /* SRC_INCLUDE_PHYSICSDATA_HPP */
