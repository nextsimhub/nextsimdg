/*!
 * @file DummyIceThermodynamics.hpp
 *
 * @date 24 Sep 2024
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef DUMMYICETHERMODYNAMICS_HPP
#define DUMMYICETHERMODYNAMICS_HPP

#include "include/IIceThermodynamics.hpp"

namespace Nextsim {

class DummyIceThermodynamics : public IIceThermodynamics {
public:
    DummyIceThermodynamics()
        : IIceThermodynamics()
    {
    }
    ~DummyIceThermodynamics() = default;

    void setData(const ModelState::DataMap& ms) override { IIceThermodynamics::setData(ms); }
    void update(const TimestepTime& tsTime) override { }
};

} /* namespace Nextsim */

#endif /* DUMMYICETHERMODYNAMICS_HPP */
