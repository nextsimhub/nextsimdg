/*!
 * @file ScopedTimer.hpp
 *
 * @date Oct 28, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef SRC_INCLUDE_SCOPEDTIMER_HPP
#define SRC_INCLUDE_SCOPEDTIMER_HPP

#include "Chrono.hpp"
#include "Timer.hpp"

#include <iostream>
#include <string>

namespace Nextsim {

class ScopedTimer {
public:
    ScopedTimer();
    ScopedTimer(const std::string& name);
    ~ScopedTimer();

    static void setTimerAddress(Timer*);
    void substitute(const std::string& newName);
    static Timer& timer();

private:
    static Timer* p_timer;
};

} /* namespace Nextsim */

#endif /* SRC_INCLUDE_SCOPEDTIMER_HPP */
