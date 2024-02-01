/*!
 * @file StressUpdateStep.hpp
 *
 * @date Feb 1, 2024
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef STRESSUPDATESTEP_HPP
#define STRESSUPDATESTEP_HPP

#include "DynamicsParameters.hpp"

namespace Nextsim {

class StressUpdateStep {
public:
    StressUpdateStep() = default;
    virtual ~StressUpdateStep() = default;
    virtual void stressUpdateHighOrder(const DynamicsParameters&) = 0;
};

} /* namespace Nextsim */

#endif /* STRESSUPDATESTEP_HPP */
