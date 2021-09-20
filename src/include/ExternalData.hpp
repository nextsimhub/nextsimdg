/*!
 * @file ExternalData.hpp
 * @date Sep 13, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef SRC_INCLUDE_EXTERNALDATA_HPP
#define SRC_INCLUDE_EXTERNALDATA_HPP

#include "BaseElementData.hpp"

namespace Nextsim {

class ExternalData: public BaseElementData {
public:
    ExternalData();
    virtual ~ExternalData();

    //! Air temperature at 2 m [˚C]
    inline double& airTemperature()
    {
        return m_tair;
    };
    inline const double& airTemperature() const
    {
        return m_tair;
    }

    //! Sea level atmospheric pressure [Pa]
    inline double& airPressure()
    {
        return m_slp;
    };
    inline const double& airPressure() const
    {
        return m_slp;
    }

    //! Water vapour mixing ratio [kg kg⁻¹]
    inline double& mixingRatio()
    {
        return m_mixrat;
    };
    inline const double& mixingRatio() const
    {
        return m_mixrat;
    }

    //! Does the element have a valid value of water vapour mixing ratio?
    inline bool hasMixingRatio() const
    {
        return (m_mixrat >= 0) && (m_mixrat <= 1);
    };

    inline double& incomingShortwave()
    {
        return m_Qsw_in;
    }
    inline const double& incomingShortwave() const
    {
        return m_Qsw_in;
    }

    inline double& incomingLongwave()
    {
        return m_Qlw_in;
    }
    inline const double& incomingLongwave() const
    {
        return m_Qlw_in;
    }
private:
    double m_tair;
    double m_slp;
    double m_mixrat;
    double m_Qsw_in;
    double m_Qlw_in;
};

} /* namespace Nextsim */

#endif /* SRC_INCLUDE_EXTERNALDATA_HPP */
