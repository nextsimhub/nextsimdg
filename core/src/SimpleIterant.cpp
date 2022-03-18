/*!
 * @file SimpleIterant.cpp
 * @date 12 Aug 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/SimpleIterant.hpp"

#include <ctime>
#include <iostream>
#include <sstream>
#include <string>

namespace Nextsim {

SimpleIterant::SimpleIterant()
{
    // It's so simple, there's nothing here
}

void SimpleIterant::init() { std::cout << "SimpleIterant::init" << std::endl; }

void SimpleIterant::start(const TimePoint& startTime)
{
    std::cout << "SimpleIterant::start at " << stringFromTimePoint(startTime) << std::endl;
}

void SimpleIterant::iterate(const TimestepTime& tst)
{
    std::cout << "SimpleIterant::iterate for " << count(tst.step) << std::endl;
}

void SimpleIterant::stop(const TimePoint& stopTime)
{
    std::cout << "SimpleIterant::stop at " << stringFromTimePoint(stopTime) << std::endl;
}

std::string SimpleIterant::stringFromTimePoint(
    const std::chrono::time_point<std::chrono::system_clock>& t)
{
    std::time_t t_c = Nextsim::Iterator::Clock::to_time_t(t);
    return std::string(ctime(&t_c));
}

std::string SimpleIterant::stringFromTimePoint(const int t)
{
    std::stringstream ss;
    ss << t;
    return ss.str();
}

template <>
std::chrono::time_point<Iterator::Clock>
SimpleIterant::zeroTime<std::chrono::time_point<Iterator::Clock>>()
{
    return Iterator::Clock::now();
}

template <> int SimpleIterant::zeroTime<int>() { return 0; }
} /* namespace Nextsim */
