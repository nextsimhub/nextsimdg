/*!
 * @file PrognosticGenerator.h
 *
 * @date Jan 28, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef CORE_SRC_PROGNOSTICGENERATOR_HPP_
#define CORE_SRC_PROGNOSTICGENERATOR_HPP_

#include "include/IPrognosticUpdater.hpp"

#include <vector>

namespace Nextsim {

class PrognosticGenerator : public IPrognosticUpdater {
public:
    PrognosticGenerator()
        : PrognosticGenerator(1)
    {
    }

    PrognosticGenerator(int nLayers)
        : m_hice(0)
        , m_cice(0)
        , m_hsnow(0)
        , m_tice(nLayers, 0.)
        , m_sst(0.)
        , m_sss(0.)
    {
    };

    virtual ~PrognosticGenerator() = default;

    PrognosticGenerator& hice(double hice)
    {
        m_hice = hice;
        return *this;
    };

    PrognosticGenerator& cice(double cice)
    {
        m_cice = cice;
        return *this;
    };

    PrognosticGenerator& hsnow(double hsnow)
    {
        m_hsnow = hsnow;
        return *this;
    };

    PrognosticGenerator& tice(const std::vector<double>& tice)
    {
        m_tice = tice;
        return *this;
    }

    PrognosticGenerator& sst(double sst)
    {
        m_sst = sst;
        return *this;
    }

    PrognosticGenerator& sss(double sss)
    {
        m_sss = sss;
        return *this;
    }

    double updatedIceThickness() const override { return m_hice; }

    double updatedSnowThickness() const override { return m_hsnow; }

    double updatedIceConcentration() const override { return m_cice; }

    const std::vector<double>& updatedIceTemperatures() const override { return m_tice; }

    double seaSurfaceTemperature() const { return m_sst; };

    double seaSurfaceSalinity() const { return m_sss; };

private:
    double m_hice;
    double m_cice;
    double m_hsnow;
    std::vector<double> m_tice;

    double m_sst;
    double m_sss;
};

} /* namespace Nextsim */

#endif /* CORE_SRC_PROGNOSTICGENERATOR_HPP_ */
