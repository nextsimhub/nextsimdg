/*!
 * @file PrognosticData.hpp
 *
 * @date Mar 1, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef CORE_SRC_INCLUDE_PROGNOSTICDATA_HPP
#define CORE_SRC_INCLUDE_PROGNOSTICDATA_HPP

#include "include/ModelModule.hpp"

#include "include/Configured.hpp"
#include "include/IceGrowth.hpp"
#include "include/Time.hpp"

namespace Nextsim {

class PrognosticData : public ModelModule, public Configured<PrognosticData> {
public:
    PrognosticData();
    virtual ~PrognosticData() = default;

    void configure() override;

    std::string getName() const override { return "PrognosticData"; };

    void setData(const ModelState& ms) override;
    ModelState getState() const override { return ModelState(); }
    ModelState getState(const OutputLevel& lvl) const override { return getState(); }

    std::set<std::string> hFields() const override { return { "h_ice", "c_cice", "h_snow" }; };
    std::set<std::string> uFields() const override { return { "u" }; }
    std::set<std::string> vFields() const override { return { "v" }; }
    std::set<std::string> zFields() const override { return { "tice" }; }

    /*!
     *  @brief Updates the state of the prognostic data for this timestep
     *
     *  @param tsInitialTime the time at the start of the timestep
     */
    void update(const TimestepTime& tsTime);

    const HField& iceThickness() { return m_thick; }

    const HField& iceConcentration() { return m_conc; }

    const HField& snowThickness() { return m_snow; }

    const UField& u() { return m_u; }

    const VField& v() { return m_v; }

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
};

} /* namespace Nextsim */

#endif /* CORE_SRC_INCLUDE_PROGNOSTICDATA_HPP */
