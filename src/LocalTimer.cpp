/*
 * @file LocalTimer.cpp
 *
 * @date Oct 28, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/LocalTimer.hpp"

namespace Nextsim {
Timer* LocalTimer::timer;

LocalTimer::LocalTimer()
    : LocalTimer("")
{
}

LocalTimer::LocalTimer(const std::string& name)
{
    timer->tick(name);
}

LocalTimer::~LocalTimer()
{
    timer->tock();
}

void LocalTimer::setTimer(Timer* timer)
{
    LocalTimer::timer = timer;
}
}
