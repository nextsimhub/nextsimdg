/*!
 * @file SimpleOutput.cpp
 *
 * @date May 25, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/SimpleOutput.hpp"
#include "include/StructureFactory.hpp"

#include <iostream>

namespace Nextsim {

void SimpleOutput::outputState(const ModelState& state, const TimestepTime& tst) const
{
    std::string timeFileName = std::to_string(tst.start) + ".nc";
    std::cout << "Outputting " << state.size() << " fields to " << timeFileName << std::endl;

    StructureFactory::fileFromState(state, timeFileName);
}
} /* namespace Nextsim */
