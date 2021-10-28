/*
 * @file LocalTimer.cpp
 *
 * @date Oct 28, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/LocalTimer.hpp"

namespace Nextsim {
Timer* LocalTimer::p_timer;

LocalTimer::LocalTimer()
    : LocalTimer("")
{
}

LocalTimer::LocalTimer(const std::string& name)
{
    p_timer->tick(name);
}

LocalTimer::~LocalTimer()
{
    p_timer->tock();
}

void LocalTimer::setTimerAddress(Timer* timer)
{
    p_timer = timer;
}

Timer& LocalTimer::timer()
{
    return *p_timer;
}
}
