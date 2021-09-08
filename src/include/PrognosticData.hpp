/*!
 * @file PrognosticData.hpp
 * @date Sep 8, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef SRC_INCLUDE_PROGNOSTICDATA_HPP
#define SRC_INCLUDE_PROGNOSTICDATA_HPP

#include <BaseElementData.hpp>

namespace Nextsim {

class PrognosticData: public BaseElementData {
public:
    PrognosticData( );
    virtual ~PrognosticData( );

    const double& iceThickness( ); //!< Effective Ice thickness [m]
    const double& iceConcentration( ); //!< Ice concentration [1]
    const double& seaSurfaceTemperature( ); //!< Sea surface temperature [˚C]
    const double& seaSurfaceSalinity( ); //!< Sea surface salinity [psu]

private:
    double m_thick; //!< Effective Ice thickness [m]
    double m_conc; //!< Ice concentration [1]
    double m_sst; //!< Sea surface temperature [˚C]
    double m_sss; //!< Sea surface salinity [psu]
};

} /* namespace Nextsim */

#endif /* SRC_INCLUDE_PROGNOSTICDATA_HPP */
