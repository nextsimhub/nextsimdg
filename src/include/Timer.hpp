/*!
 * @file Timer.hpp
 * @date 11 Aug 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef SRC_INCLUDE_TIMER_HPP
#define SRC_INCLUDE_TIMER_HPP

#include <string>

namespace Nextsim {

class Timer {
public:
    Timer();

    void tick(const std::string& timerName);
    void tock(const std::string& timerName);
    std::string report() const;
};

} /* namespace Nextsim */

#endif /* SRC_INCLUDE_TIMER_HPP */
