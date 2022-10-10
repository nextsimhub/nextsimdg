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

//! A class implementing a constant value AtmosphereState.
class DummyAtmosphereState : public AtmosphereState {
public:
    DummyAtmosphereState() = default;
    ~DummyAtmosphereState() = default;

    void setData(const ModelState::DataMap&) override;
    std::string getName() const override { return "DummyAtmosphereState"; }

protected:
    //! Performs the implementation specific updates. Does nothing.
    void updateSpecial(const TimestepTime&) override { }
};

} /* namespace Nextsim */

#endif /* DUMMYATMOSPHERESTATE_HPP */
