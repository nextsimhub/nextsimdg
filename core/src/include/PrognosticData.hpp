/*!
 * @file PrognosticData.hpp
 *
 * @date Mar 1, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef PROGNOSTICDATA_HPP
#define PROGNOSTICDATA_HPP

#include "ModelComponent.hpp"
#include "include/Configured.hpp"
#include "include/IAtmosphereBoundary.hpp"
#include "include/IceGrowth.hpp"
#include "include/Time.hpp"

namespace Nextsim {

/*!
 * A class defining a configurable ModelComponent that stores the prognostic
 * data values and handles their updates in the timestep, including all calls
 * to the variables those calculations depend on.
 */
class PrognosticData : public ModelComponent, public Configured<PrognosticData> {
public:
    PrognosticData();
    virtual ~PrognosticData() = default;

    void configure() override;

    std::string getName() const override { return "PrognosticData"; };

    void setData(const ModelState::DataMap& ms) override;
    ModelState getState() const override;
    ModelState getState(const OutputLevel& lvl) const override { return getState(); }
    ModelState getStateRecursive(const OutputSpec& os) const override;

    static HelpMap& getHelpText(HelpMap& map, bool getAll);
    static HelpMap& getHelpRecursive(HelpMap& map, bool getAll);

    std::unordered_set<std::string> hFields() const override
    {
        return { "h_ice", "c_cice", "h_snow" };
    };
    std::unordered_set<std::string> uFields() const override { return { "u" }; }
    std::unordered_set<std::string> vFields() const override { return { "v" }; }
    std::unordered_set<std::string> zFields() const override { return { "tice" }; }

    /*!
     *  @brief Updates the state of the prognostic data for this timestep
     *
     *  @param tsInitialTime the time at the start of the timestep
     */
    void update(const TimestepTime& tsTime);

    //! Returns a const reference to the cell-averaged ice thickness field.
    const HField& iceThickness() { return m_thick; }

    //! Returns a reference to the ice concentration field
    const HField& iceConcentration() { return m_conc; }

    //! Returns a const reference to the cell-averaged snow thickness field.
    const HField& snowThickness() { return m_snow; }

    //! Returns a const reference to the eastward component of the ice drift velocity.
    const UField& u() { return m_u; }

    //! Returns a const reference to the northward component of the ice drift velocity.
    const VField& v() { return m_v; }

    //! Returns a const reference to the (three dimensional) ice temperature field.
    const ZField& iceTemperature() { return m_tice; }

private:
    HField m_thick;
    HField m_conc;
    ZField m_tice;
    HField m_snow;
    UField m_u;
    VField m_v;
    double m_dt;

    IceGrowth iceGrowth;
    IAtmosphereBoundary *pAtmBdy;
};

} /* namespace Nextsim */

#endif /* PROGNOSTICDATA_HPP */
