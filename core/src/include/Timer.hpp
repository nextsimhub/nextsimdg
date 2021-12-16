/*!
 * @file Timer.hpp
 * @date 11 Aug 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef SRC_INCLUDE_TIMER_HPP
#define SRC_INCLUDE_TIMER_HPP

#include "Chrono.hpp"

#include <chrono>
#include <ctime>
#include <forward_list>
#include <map>
#include <ostream>
#include <set>
#include <string>

namespace Nextsim {

//! A class for a hierarchical timer functions.
class Timer {
public:
    typedef std::string Key;

    typedef Chrono::WallTimePoint WallTimePoint;
    typedef Chrono::WallTimeDuration WallTimeDuration;
    typedef Chrono::CpuTimePoint CpuTimePoint;
    typedef Chrono::CpuTimeDuration CpuTimeDuration;

    typedef std::forward_list<Key> TimerPath;

    //! Creates a Timer with an unnamed root node.
    Timer();
    /*!
     * @brief Creates a Timer with a named root node.
     *
     * @param rootKey Name of the root node.
     */
    Timer(const Key& rootKey);
    virtual ~Timer() = default;

    /*!
     * @brief Starts a named timer.
     *
     * @param timerName Name of the timer to be started.
     */
    void tick(const Key& timerName);
    /*!
     * @brief Stops a named timer.
     *
     * @param timerName Name of the timer to be stopped.
     */
    void tock(const Key& timerName);
    //! @brief Stop the last timer to be started.
    void tock();

    /*!
     * @brief Returns the elapsed time without stopping the timer.
     *
     * @param timerName the name of the timer to interrogate.
     */
    double lap(const Key& timerName) const;
    /*!
     * @brief Returns the elapsed time.
     *
     * @param timerName the name of the timer to interrogate.
     */
    double elapsed(const Key& timerName) const;

    /*!
     * @brief Prints the status of a named timer to an ostream.
     *
     * @param timerName The timer to be printed.
     * @param os The ostream to print to.
     */
    std::ostream& report(const Key& timerName, std::ostream& os) const;
    /*!
     * @brief Prints the status of all the timers to an ostream.
     *
     * @param os The ostream to print to.
     */
    std::ostream& report(std::ostream& os) const;
    /*!
     * @brief Prints the status of a timer specified by a path to an ostream.
     *
     * @param path The path to the timer to be printed.
     * @param os The ostream to print to.
     */
    std::ostream& report(const TimerPath&, std::ostream& os) const;

    /*!
     * @brief Adds an additional time increment to a timer.
     *
     * @param path Path to the timer to be incremented.
     * @param additionalWall Wall clock duration to be added.
     * @param additionalCpu CPU clock duration to be added.
     * @param additionalTicks Activation ticks to be added.
     */
    void additionalTime(const TimerPath& path, WallTimeDuration additionalWall,
        CpuTimeDuration additionalCpu, int additionalTicks);
    //! Returns the timer path to the currently running timer.
    TimerPath currentTimerNodePath() const;

    //! Deletes all timers except the root, which is reset.
    void reset();
    //! Static timer for general use.
    static Timer main;

private:
    struct TimerNode {
        TimerNode();
        Key name;

        Chrono timeKeeper;

        std::map<Key, TimerNode> childNodes;
        TimerNode* parent;

        void tick();
        void tock();
        std::ostream& report(std::ostream& os, const std::string& prefix) const;
        std::ostream& reportAll(std::ostream& os, const std::string& prefix) const;

        TimerPath searchDescendants(const Key& timerName) const;
    };

    TimerPath pathToFirstMatch(const Key&) const;

    TimerNode root;
    TimerNode* current;
};

} /* namespace Nextsim */

//! Overload of the output operator for a Timer.
std::ostream& operator<<(std::ostream& os, const Nextsim::Timer& tim);

#endif /* SRC_INCLUDE_TIMER_HPP */
