/*!
 * @author  Tim Spain <timothy.spain@nersc.no>
 */

#include "include/Iterator.hpp"

#include <sstream>

namespace Nextsim {

void Iterator::setStartStopStep(TimePoint startTime, TimePoint stopTime, Duration timestep)
{
    this->startTime = startTime;
    this->stopTime = stopTime;
    this->timestep = timestep;
}

void Iterator::run()
{
    iterant.start(startTime);

    for (auto t = startTime; t < stopTime; t += timestep) {
        TimestepTime tsTime = { t, timestep };
        try {
            iterant.iterate(tsTime);
        } catch (const std::exception& e) {
            iterant.stop(t);
            throw std::runtime_error(
                e.what() + std::string(" Execution halted at time step ") + tsTime.start.format());
        }
    }

    iterant.stop(stopTime);
}

} /* namespace Nextsim */
