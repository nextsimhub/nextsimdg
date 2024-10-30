/*!
 * @file DummyDynamics.hpp
 *
 * @date 24 Sep 2024
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef DUMMYDYNAMICS_HPP
#define DUMMYDYNAMICS_HPP

#include "include/IDynamics.hpp"

#include "include/ModelArray.hpp"
#include "include/ModelComponent.hpp"

namespace Nextsim {
class DummyDynamics : public IDynamics {
public:
    DummyDynamics() = default;

    std::string getName() const override { return "DummyDynamics"; }
    void update(const TimestepTime& tst) override
    {
        uice = 0.;
        vice = 0.;
    };
protected:
    ModelArray& advectHField(ModelArray& field, const std::string& fieldName) override
    {
        return field;
    }
};
}

#endif /* DUMMYDYNAMICS_HPP */
