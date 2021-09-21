/*!
 * @file PrognosticData.hpp
 * @date Sep 8, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef SRC_INCLUDE_PROGNOSTICDATA_HPP
#define SRC_INCLUDE_PROGNOSTICDATA_HPP

#include "BaseElementData.hpp"

namespace Nextsim {

class PrognosticData: public BaseElementData {
public:
    PrognosticData( ) = default;
    ~PrognosticData( ) = default;

    //! Effective Ice thickness [m]
    inline const double& iceThickness( ) const
    {
        return m_thick;
    }

    //! Ice concentration [1]
    inline const double& iceConcentration( ) const
    {
        return m_conc;
    }

    //! Sea surface temperature [˚C]
    inline const double& seaSurfaceTemperature( ) const
    {
        return m_sst;
    }

    //! Sea surface salinity [psu]
    const double& seaSurfaceSalinity( ) const
    {
        return m_sss;
    }

private:
    double m_thick; //!< Effective Ice thickness [m]
    double m_conc; //!< Ice concentration [1]
    double m_sst; //!< Sea surface temperature [˚C]
    double m_sss; //!< Sea surface salinity [psu]
};

} /* namespace Nextsim */

#endif /* SRC_INCLUDE_PROGNOSTICDATA_HPP */
