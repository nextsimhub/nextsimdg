/*!
 * @file PDTestDynamics.hpp
 *
 * @date Jan 21, 2025
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef PDTESTDYNAMICS_HPP
#define PDTESTDYNAMICS_HPP

#include "include/DummyDynamics.hpp"

namespace Nextsim {
class PDTestDynamics : public DummyDynamics {
public:
    std::string getName() const override { return "PDTestDynamics"; }

    void setData (const ModelState::DataMap& ms) override
    {
        dataMap = ms;
    }
    ModelState getState(const OutputLevel&) const override { return getState(); }
    ModelState getState() const override
    {
        return { dataMap, {} };
    }

    ModelState getStateRecursive(const OutputSpec& os) const override
    {
        return { dataMap, {} };
    }
private:
    ModelState::DataMap dataMap;
};
}

#endif /* PDTESTDYNAMICS_HPP */
