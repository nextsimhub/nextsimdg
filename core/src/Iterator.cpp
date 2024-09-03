/*!
 * @file Iterator.cpp
 * @date 11 Aug 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/Iterator.hpp"

#include <sstream>

namespace Nextsim {

Iterator::Iterator()
{
    iterant = new NullIterant();
}

Iterator::Iterator(Iterant* iterant)
    : iterant(iterant)
{
}

void Iterator::setIterant(Iterant* iterant) { this->iterant = iterant; }

void Iterator::setStartStopStep(TimePoint startTime, TimePoint stopTime, Duration timestep)
{
    this->startTime = startTime;
    this->stopTime = stopTime;
    this->timestep = timestep;
}

TimePoint Iterator::parseAndSet(const std::string& startTimeStr, const std::string& stopTimeStr,
    const std::string& durationStr, const std::string& stepStr)
{
    std::stringstream ss(startTimeStr);
    ss >> startTime;
    ss = std::stringstream(stepStr);
    ss >> timestep;
    if (!durationStr.empty()) {
        ss = std::stringstream(durationStr);
        Duration duration;
        ss >> duration;
        stopTime = startTime + duration;
    } else {
        ss = std::stringstream(stopTimeStr);
        ss >> stopTime;
    }

    return startTime;
}

void Iterator::run()
{
    iterant->start(startTime);

    for (auto t = startTime; t < stopTime; t += timestep) {
        TimestepTime tsTime = { t, timestep };
        iterant->iterate(tsTime);
    }

    iterant->stop(stopTime);
}

} /* namespace Nextsim */
