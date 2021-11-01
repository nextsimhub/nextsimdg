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

class PhysicsData: public BaseElementData {
public:
    PhysicsData() = default;
    ~PhysicsData() = default;

    //! Rate of evaporation (+ve) or condensation (-ve) [kg s⁻¹ m⁻²]
    inline double& evaporationRate()
    {
        return m_evap;
    };
    //! Rate of sublimation ice to vapour [kg s⁻¹ m⁻²]
    inline double& sublimationRate()
    {
        return m_subl;
    }
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
    //! Specific humidity over the ice [kg kg⁻¹]
    inline double& specificHumidityIce()
    {
        return m_sphumi;
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
    //! Specific heat capacity of wet air
    inline double& heatCapacityWetAir()
    {
        return m_cspec;
    }
    //! Pressure due to wind drag [Pa]
    inline double& dragPressure()
    {
        return m_tau;
    }
    //! Total heat flux over open water [W m⁻²]
    inline double& QOpenWater()
    {
        return m_Qow;
    }
    //! Ocean to atmosphere longwave heat flux over open water [W m⁻²]
    inline double& QLongwaveOpenWater()
    {
        return m_Qlwow;
    }
    //! Ocean to atmosphere longwave heat flux over open water [W m⁻²]
    inline double& QShortwaveOpenWater()
    {
        return m_Qswow;
    }
    //! Ocean to atmosphere longwave heat flux over open water [W m⁻²]
    inline double& QLatentHeatOpenWater()
    {
        return m_Qlhow;
    }
    //! Ocean to atmosphere longwave heat flux over open water [W m⁻²]
    inline double& QSensibleHeatOpenWater()
    {
        return m_Qshow;
    }

    //! Total heat flux over ice [W m⁻²]
    inline double& QIce()
    {
        return m_Qi;
    }
    //! Ocean to atmosphere longwave heat flux over ice [W m⁻²]
    inline double& QLongwaveIce()
    {
        return m_Qlwi;
    }
    //! Ocean to atmosphere longwave heat flux over ice [W m⁻²]
    inline double& QShortwaveIce()
    {
        return m_Qswi;
    }
    //! Ocean to atmosphere longwave heat flux over ice [W m⁻²]
    inline double& QLatentHeatIce()
    {
        return m_Qlhi;
    }
    //! Ocean to atmosphere longwave heat flux over ice [W m⁻²]
    inline double& QSensibleHeatIce()
    {
        return m_Qshi;
    }
    //! Derivative of ice-atmosphere heat flux with respect to ice surface temperature [W m⁻² K⁻¹]
    inline double& QDerivativeWRTTemperature()
    {
        return m_dQ_dT;
    }

    //! Ice to ocean heat flux [W m⁻²]
    inline double& QIceOceanHeat()
    {
        return m_Qio;
    }

    //! Total ice-atmosphere heat flux [W m⁻²]
    inline double& QIceAtmosphere()
    {
        return m_Qia;
    }
    //! Mean thickness of ice (averaged over ice covered fraction) [m]
    inline double& iceTrueThickness()
    {
        return m_hi;
    }

    //! Mean thickness of snow (averaged over ice covered fraction) [m]
    inline double& snowTrueThickness()
    {
        return m_hs;
    }

    //! Updated value of the ice surface temperature [˚C]
    inline double& updatedIceSurfaceTemperature()
    {
        return m_TiceNew[0];
    }

    //! Total amount of ice that resulted from the flooding of snow [m]
    inline double& totalIceFromSnow()
    {
        return m_hifroms;
    }

    inline double& oceanAlbedo()
    {
        return m_oceanAlbedo;
    }
private:
    double m_evap;
    double m_subl;
    double m_rho;
    double m_wspeed;
    double m_sphumw;
    double m_sphumi;
    double m_sphuma;
    double m_cspec;
    double m_tau;

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

    // thermodynamic values
    double m_hi;
    double m_hs;
    std::array<double, N_ICE_TEMPERATURES> m_TiceNew;
    double m_hifroms;

    // Ocean albedo
    static double m_oceanAlbedo;
};

} /* namespace Nextsim */

#endif /* SRC_INCLUDE_PHYSICSDATA_HPP */
