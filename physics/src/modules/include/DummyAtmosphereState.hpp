/*!
 * @file DummyAtmosphereState.hpp
 *
 * @date May 10, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef DUMMYATMOSPHERESTATE_HPP
#define DUMMYATMOSPHERESTATE_HPP

#include "AtmosphereState.hpp"

namespace Nextsim {

class DummyAtmosphereState : public AtmosphereState {
public:
    DummyAtmosphereState() = default;
    ~DummyAtmosphereState() = default;

    void setData(const ModelState&) override;
    std::string getName() const override { return "DummyAtmosphereState"; }

protected:
    void updateSpecial(const TimestepTime&) override { }
};

} /* namespace Nextsim */

#endif /* DUMMYATMOSPHERESTATE_HPP */
