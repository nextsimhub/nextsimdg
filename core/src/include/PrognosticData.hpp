/*!
 * @file PrognosticData.hpp
 * @date Sep 8, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef SRC_INCLUDE_PROGNOSTICDATA_HPP
#define SRC_INCLUDE_PROGNOSTICDATA_HPP

#include "BaseElementData.hpp"
#include "Configured.hpp"
#include "include/IFreezingPoint.hpp"
#include <array>

namespace Nextsim {

const int N_ICE_TEMPERATURES = 3;

//! A class holding all of the data for an element that is carried from one
//! timestep to another.
class PrognosticData : public BaseElementData, public Configured<PrognosticData> {
public:
    PrognosticData();
    ~PrognosticData() = default;

    void configure() override;

    //! Effective Ice thickness [m]
    inline double iceThickness() const { return m_thick; }
    //! True ice thickness [m]. Zero concentration means no ice thickness
    inline double iceTrueThickness() const { return (m_conc != 0) ? m_thick / m_conc : 0; }

    //! Ice concentration [1]
    inline double iceConcentration() const { return m_conc; }

    //! Sea surface temperature [˚C]
    inline double seaSurfaceTemperature() const { return m_sst; }

    //! Sea surface salinity [psu]
    inline double seaSurfaceSalinity() const { return m_sss; }

    //! Ice temperatures [˚C]
    inline const std::array<double, N_ICE_TEMPERATURES>& iceTemperatures() const { return m_tice; }
    template <int I> double iceTemperature() const
    {
        static_assert(
            I < N_ICE_TEMPERATURES, "Ice layer indices must be 0 <= I < N_ICE_TEMPERATURES.");
        return m_tice[I];
    }

    //! Mean snow thickness [m]
    inline double snowThickness() const { return m_snow; }
    //! Mean snow thickness over ice [m]
    inline double snowTrueThickness() const { return (m_conc != 0) ? m_snow / m_conc : 0; }

    //! Salinity dependent freezing point [˚C]
    inline double freezingPoint() const { return (*m_freezer)(m_sss); }

    //! Timestep [s]
    inline double timestep() const { return m_dt; }
    //! Set a new value for the timestep
    static void setTimestep(double newDt) { m_dt = newDt; }

    /*!
     * @brief Set the data from passed arguments.
     *
     * @param h Ice thickness [m]
     * @param c Ice concentration [1]
     * @param t Sea surface temperature [˚C]
     * @param s Sea surface salinity [PSU]
     * @param hs Snow thickness [m]
     * @param tice Array of ice temperatures [˚C]
     */
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
