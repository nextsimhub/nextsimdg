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
    ModelState getStateDiagnostic() const override
    {
        return { dataMap, {} };
    }

    ModelState getStatePrognostic() const override
    {
        return { dataMap, {} };
    }
private:
    ModelState::DataMap dataMap;
};
}

#endif /* PDTESTDYNAMICS_HPP */
