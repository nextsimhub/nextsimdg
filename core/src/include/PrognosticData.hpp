/*!
 * @file PrognosticData.hpp
 *
 * @date Mar 1, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef CORE_SRC_INCLUDE_PROGNOSTICDATA_HPP
#define CORE_SRC_INCLUDE_PROGNOSTICDATA_HPP

#include "include/ModelModule.hpp"

namespace Nextsim {

class PrognosticData : public ModelModule {
public:
    PrognosticData();
    virtual ~PrognosticData() = default;

    std::string getName() const override { return "PrognosticData"; };

    void setData(const ModelState& ms) override { }
    ModelState getState() const override { return ModelState(); }
    ModelState getState(const OutputLevel& lvl) const override { return getState(); }

    HField& iceThickness() { return m_thick; }
    double& iceThickness(size_t i, size_t j) { return m_thick(i, j); }

    HField& iceConcentration() { return m_conc; }
    double& iceConcentration(size_t i, size_t j) { return m_conc(i, j); }

    HField& snowThickness() { return m_snow; }
    double& snowThickness(size_t i, size_t j) { return m_snow(i, j); }

    UField& u() { return m_u; }
    double& u(size_t i, size_t j) { return m_u(i, j); }

    VField& v() { return m_v; }
    double& v(size_t i, size_t j) { return m_v(i, j); }

    ZField& iceTemperature() { return m_tice; }

private:
    HField m_thick;
    HField m_conc;
    ZField m_tice;
    HField m_snow;
    UField m_u;
    VField m_v;
    double m_dt;
};

} /* namespace Nextsim */

#endif /* CORE_SRC_INCLUDE_PROGNOSTICDATA_HPP */
