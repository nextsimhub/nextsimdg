/*!
 * @file FiniteElementAtmos.hpp
 *
 * @date May 2, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef FINITEELEMENTATMOS_HPP
#define FINITEELEMENTATMOS_HPP

#include "IAtmosphericState.hpp"

namespace Nextsim {

class FiniteElementAtmos: public IAtmosphericState {
public:
    FiniteElementAtmos()
        : IAtmosphericState()
    {

    }

    std::string getName() const override { return "FiniteElementAtmos"; }

    void update(const TimestepTime&) override;
    void updateElement(size_t i, const TimestepTime& tst);
};

} /* namespace Nextsim */

#endif /* FINITEELEMENTATMOS_HPP */
