/*
 * @file LocalTimer.hpp
 *
 * @date Oct 28, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef SRC_INCLUDE_LOCALTIMER_HPP
#define SRC_INCLUDE_LOCALTIMER_HPP

#include "Chrono.hpp"
#include "Timer.hpp"

#include <iostream>
#include <string>

namespace Nextsim {

class LocalTimer {
public:
    LocalTimer();
    LocalTimer(const std::string& name);
    ~LocalTimer();

    static void setTimer(Timer* timer);
private:
    static Timer *timer;
};

} /* namespace Nextsim */

#endif /* SRC_INCLUDE_LOCALTIMER_HPP */
