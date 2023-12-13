/*!
 * @file DummyIceThermodynamics.hpp
 *
 * @date 18 Apr 2023
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef DUMMYICETHERMODYNAMICS_HPP
#define DUMMYICETHERMODYNAMICS_HPP

#include "include/IIceThermodynamics.hpp"

#include "include/NZLevels.hpp"

#include "include/NZLevels.hpp"

namespace Nextsim {

class DummyIceThermodynamics : public IIceThermodynamics {
public:
    DummyIceThermodynamics()
    : IIceThermodynamics()
    {
        NZLevels::set(getNZLevels());
    }
    ~DummyIceThermodynamics() = default;

    ModelState getStateRecursive(const OutputSpec& os) const override
    {
        return ModelState();
    }

    void setData(const ModelState::DataMap& ms) override
    {
        IIceThermodynamics::setData(ms);
    }
    void update(const TimestepTime& tsTime) override {}

    size_t getNZLevels() const override
    {
        return 1;
    } // 1 is the minimum, I guess
};

} /* namespace Nextsim */

#endif /* DUMMYICETHERMODYNAMICS_HPP */
