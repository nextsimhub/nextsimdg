/*!
 * @file DummyIceSpread.hpp
 *
 * @date 04 Jun 2025
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef DUMMYICESPREAD_HPP
#define DUMMYICESPREAD_HPP

#include "include/ILateralIceSpread.hpp"

namespace Nextsim {

class DummyIceSpread : public ILateralIceSpread {
public:
    DummyIceSpread()
        : ILateralIceSpread()
    {
    }
    ~DummyIceSpread() = default;

    void update(const TimestepTime& tstep) override {
        // No-op for dummy implementation
    };
};

} /* namespace Nextsim */

#endif /* DUMMYICESPREAD_HPP */
