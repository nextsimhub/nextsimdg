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

//! A class providing a timer aware of the calling context
class ScopedTimer {
public:
    ScopedTimer();
    //! Creates a scoped timer with a name
    ScopedTimer(const std::string& name);
    ~ScopedTimer();

    /*!
     * Sets the address of the timer which provides the timing functions
     *
     * @param timer A pointer to an instance of the Timer class.
     */
    static void setTimerAddress(Timer* timer);
    /*!
     * Substitutes the currently running timer with a new one with a different
     * name.
     *
     * @param newName the name of the new timer
     */
    void substitute(const std::string& newName);
    //! Returns a reference to the underlying Timer.
    static Timer& timer();

private:
    static Timer* p_timer;
};

} /* namespace Nextsim */

#endif /* SRC_INCLUDE_SCOPEDTIMER_HPP */
