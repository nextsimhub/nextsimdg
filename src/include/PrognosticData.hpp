/*!
 * @file PrognosticData.hpp
 * @date Sep 8, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef SRC_INCLUDE_PROGNOSTICDATA_HPP
#define SRC_INCLUDE_PROGNOSTICDATA_HPP

#include "BaseElementData.hpp"
#include "IFreezingPoint.hpp"
#include <array>

namespace Nextsim {

const int N_ICE_TEMPERATURES = 3;

class PrognosticData : public BaseElementData {
public:
    PrognosticData();
    ~PrognosticData() = default;

    //! Effective Ice thickness [m]
    inline const double& iceThickness() const { return m_thick; }

    //! Ice concentration [1]
    inline const double& iceConcentration() const { return m_conc; }

    //! Sea surface temperature [˚C]
    inline const double& seaSurfaceTemperature() const { return m_sst; }

    //! Sea surface salinity [psu]
    inline const double& seaSurfaceSalinity() const { return m_sss; }

    //! Ice temperatures [˚C]
    inline const std::array<double, N_ICE_TEMPERATURES>& iceTemperatures() const { return m_tice; }

    //! Mean snow thickness [m]
    inline const double& snowThickness() const { return m_snow; }
    //! Salinity dependent freezing point
    inline const double freezingPoint() const { return (*m_freezer)(m_sss); }

    //! Timestep [s]
    inline const double& timestep() const { return m_dt; }
    //! Set a new value for the timestep
    static void setTimestep(double newDt) { m_dt = newDt; }

    inline static PrognosticData generate(double h, double c, double t, double s, double hs,
        std::array<double, N_ICE_TEMPERATURES> tice)
    {
        PrognosticData data;
        data.m_thick = h;
        data.m_conc = c;
        data.m_sst = t;
        data.m_sss = s;
        data.m_snow = hs;
        data.m_tice = tice;

        return data;
    }

private:
    double m_thick; //!< Effective Ice thickness [m]
    double m_conc; //!< Ice concentration [1]
    double m_sst; //!< Sea surface temperature [˚C]
    double m_sss; //!< Sea surface salinity [psu]
    std::array<double, N_ICE_TEMPERATURES> m_tice; //!< Ice temperature [˚C]
    double m_snow; //!< Mean snow thickness [m]

    static double m_dt; //!< Current timestep, shared by all elements
    static IFreezingPoint* m_freezer;
};

} /* namespace Nextsim */

#endif /* SRC_INCLUDE_PROGNOSTICDATA_HPP */
