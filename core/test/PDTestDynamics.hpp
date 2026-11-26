/*!
 * @author  Tim Spain <timothy.spain@nersc.no>
 */

#ifndef PDTESTDYNAMICS_HPP
#define PDTESTDYNAMICS_HPP

#include "include/DummyDynamics.hpp"

namespace Nextsim {
class PDTestDynamics : public DummyDynamics {
public:
    std::string getName() const override { return "PDTestDynamics"; }

    void setData(const ModelState::DataMap& ms) override { dataMap = ms; }

    void appendStatePrognostic(ModelState& state) const override { state.data = dataMap; }

private:
    ModelState::DataMap dataMap;
};
}

#endif /* PDTESTDYNAMICS_HPP */
