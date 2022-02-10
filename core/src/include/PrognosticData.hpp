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
#include "include/IPrognosticUpdater.hpp"
#include "include/PrognosticGenerator.hpp"

#include <array>
#include <vector>

namespace Nextsim {

//! A class holding all of the data for an element that is carried from one
//! timestep to another.
class PrognosticData : public BaseElementData, public Configured<PrognosticData> {
public:
    //! Constructs an instance with a default of 1 ice layer.
    PrognosticData();
    //! Constructs an instance with a number of ice layers.
    PrognosticData(int nIceLayers);
    PrognosticData(const PrognosticGenerator&);
    ~PrognosticData() = default;

    /*!
     * Assigns directly from an IPrognosticUpdater.
     * @param up the updater containing the new data values.
     */
    PrognosticData& operator=(const PrognosticGenerator& up);

    /*!
     * Assigns new values from an implementation of IPrognosticUpdater
     * @param updater The IPrognosticUpdater providing the updated values
     */
    PrognosticData& updateAndIntegrate(const IPrognosticUpdater& updater);

    /*!
     * Sets the sea surface parameters, if required
     * @param sst sea surface temperature
     * @param sss sea surface salinity
     */
    PrognosticData& setSeaSurface(double sst, double sss);

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
    inline const std::vector<double>& iceTemperatures() const { return m_tice; }
    template <int I> double iceTemperature() const { return m_tice[I]; }
    double iceTemperature(int i) const { return m_tice[i]; };

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

    //! Returns the number of ice layers in this element.
    int nIceLayers() const { return m_tice.size(); };

private:
    double m_thick; //!< Effective Ice thickness [m]
    double m_conc; //!< Ice concentration [1]
    double m_sst; //!< Sea surface temperature [˚C]
    double m_sss; //!< Sea surface salinity [psu]
    std::vector<double> m_tice; //!< Ice temperature [˚C]
    double m_snow; //!< Mean snow thickness [m]

    static double m_dt; //!< Current timestep, shared by all elements
    static IFreezingPoint* m_freezer;

    static void copyInIceLayerData(const std::vector<double>& src, std::vector<double>& tgt);
};

} /* namespace Nextsim */

#endif /* SRC_INCLUDE_PROGNOSTICDATA_HPP */
