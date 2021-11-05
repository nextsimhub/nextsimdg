/*!
 * @file ScopedTimer.cpp
 *
 * @date Oct 28, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/ScopedTimer.hpp"

namespace Nextsim {
Timer* ScopedTimer::p_timer;

ScopedTimer::ScopedTimer()
    : ScopedTimer("")
{
}

ScopedTimer::ScopedTimer(const std::string& name) { p_timer->tick(name); }

ScopedTimer::~ScopedTimer() { p_timer->tock(); }

void ScopedTimer::substitute(const std::string& newName)
{
    p_timer->tock();
    p_timer->tick(newName);
}

void ScopedTimer::setTimerAddress(Timer* timer) { p_timer = timer; }

Timer& ScopedTimer::timer() { return *p_timer; }
}
