/*!
 * @file DummyOceanState.hpp
 *
 * @date May 10, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef DUMMYOCEANSTATE_HPP
#define DUMMYOCEANSTATE_HPP

#include "include/OceanState.hpp"

namespace Nextsim {

//! A class implementing a constant value OceanState.
class DummyOceanState : public OceanState {
public:
    DummyOceanState() = default;
    ~DummyOceanState() = default;

    void setData(const ModelState::DataMap&) override;
    std::string getName() const override { return "DummyOceanState"; }

protected:
    //! Performs the implementation specific updates. Does nothing.
    void updateSpecial(const TimestepTime&) override { }
};

} /* namespace Nextsim */

#endif /* DUMMYOCEANSTATE_HPP */
