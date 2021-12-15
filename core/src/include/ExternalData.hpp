/*!
 * @file ExternalData.hpp
 * @date Sep 13, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef SRC_INCLUDE_EXTERNALDATA_HPP
#define SRC_INCLUDE_EXTERNALDATA_HPP

#include "BaseElementData.hpp"
#include "constants.hpp"

namespace Nextsim {

//! A class holding all of the data for an element that is imported from
// external sources (coupled models, climatologies).
class ExternalData : public BaseElementData {
public:
    ExternalData() = default;
    ~ExternalData() = default;

    //! Reference to the air temperature at 2 m [˚C]
    inline double& airTemperature() { return m_tair; };
    //! Air temperature at 2 m [˚C]
    inline double airTemperature() const { return m_tair; }

    //! Reference to the dew point temperature at 2 m [˚C]
    inline double& dewPoint2m() { return m_dair; };
    //! Dew point temperature at 2 m [˚C]
    inline double dewPoint2m() const { return m_dair; };

    //! Reference to the sea level atmospheric pressure [Pa]
    inline double& airPressure() { return m_slp; };
    //! Sea level atmospheric pressure [Pa]
    inline double airPressure() const { return m_slp; }

    //! Reference to the water vapour mixing ratio [kg kg⁻¹]
    inline double& mixingRatio() { return m_mixrat; };
    //! Water vapour mixing ratio [kg kg⁻¹]
    inline double mixingRatio() const { return m_mixrat; }

    //! Does the element have a valid value of water vapour mixing ratio?
    inline bool hasMixingRatio() const { return (m_mixrat >= 0) && (m_mixrat <= 1); };

    //! Reference to the incoming short wave radiation flux [W m⁻²]
    inline double& incomingShortwave() { return m_Qsw_in; }
    //! Incoming short wave radiation flux [W m⁻²]
    inline double incomingShortwave() const { return m_Qsw_in; }

    //! Reference to the incoming long wave radiation flux [W m⁻²]
    inline double& incomingLongwave() { return m_Qlw_in; }
    //! Incoming long wave radiation flux [W m⁻²]
    inline double incomingLongwave() const { return m_Qlw_in; }

    //! Reference to the depth of the ocean mixed layer [m]
    inline double& mixedLayerDepth() { return m_mld; };
    //! Depth of the ocean mixed layer [m]
    inline double mixedLayerDepth() const { return m_mld; }
    //! The areal mixed layer heat capacity [J K⁻¹ m⁻²]
    inline double mixedLayerBulkHeatCapacity() const { return m_mld * Water::rhoOcean * Water::cp; }

    //! Reference to the snowfall rate [kg m⁻² s⁻¹]
    inline double& snowfall() { return m_snowfall; }
    //! Snowfall rate [kg m⁻² s⁻¹]
    inline double snowfall() const { return m_snowfall; }

private:
    double m_tair;
    double m_dair;
    double m_slp;
    double m_mixrat;
    double m_Qsw_in;
    double m_Qlw_in;
    double m_mld;

    double m_snowfall;
};

} /* namespace Nextsim */

#endif /* SRC_INCLUDE_EXTERNALDATA_HPP */
