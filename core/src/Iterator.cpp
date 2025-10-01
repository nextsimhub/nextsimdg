/*!
 * @author  Tim Spain <timothy.spain@nersc.no>
 * @author  Joe Wallwork <jw2423@cam.ac.uk>
 */

#include "include/Iterator.hpp"

#include <sstream>

namespace Nextsim {

void Iterator::setStartStopStep(TimePoint startTime, TimePoint stopTime, Duration timestep)
{
    TimePartition& timePartition = Iterator::getTimePartition();
    timePartition.startTime = startTime;
    timePartition.stopTime = stopTime;
    timePartition.timestep = timestep;
}

TimePoint Iterator::parseAndSet(const std::string& startTimeStr, const std::string& stopTimeStr,
    const std::string& durationStr, const std::string& stepStr)
{
    TimePartition& timePartition = Iterator::getTimePartition();

    std::stringstream ss(startTimeStr);
    ss >> timePartition.startTime;
    ss = std::stringstream(stepStr);
    ss >> timePartition.timestep;
    if (!durationStr.empty()) {
        ss = std::stringstream(durationStr);
        Duration duration;
        ss >> duration;
        timePartition.stopTime = timePartition.startTime + duration;
    } else {
        ss = std::stringstream(stopTimeStr);
        ss >> timePartition.stopTime;
    }

    return timePartition.startTime;
}

void Iterator::run()
{
    TimePartition& timePartition = Iterator::getTimePartition();

    iterant.start(timePartition.startTime);

    for (auto t = timePartition.startTime; t < timePartition.stopTime;
         t += timePartition.timestep) {
        TimestepTime tsTime = { t, timePartition.timestep };
        try {
            iterant.iterate(tsTime);
        } catch (const std::exception& e) {
            iterant.stop(t);
            throw std::runtime_error(
                e.what() + std::string(" Execution halted at time step ") + tsTime.start.format());
        }
    }

    iterant.stop(timePartition.stopTime);
}

} /* namespace Nextsim */
