/*!
 * @file Timed.cpp
 * @date 11 Aug 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include <sstream>
#include <string>

#include "Timed.hpp"
#include "Timer.hpp"

namespace Nextsim {
Timer Timed::timer;

Timed::Timed() { }

void Timed::tick(const std::string& timerName) { timer.tick(timerName); }

void Timed::tock(const std::string& timerName) { timer.tock(timerName); }

std::string Timed::report()
{
    std::stringstream ss;
    timer.report(ss);
    return ss.str();
}
}
