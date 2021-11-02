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

class Timer {
public:
    typedef std::string Key;

    typedef Chrono::WallTimePoint WallTimePoint;
    typedef Chrono::WallTimeDuration WallTimeDuration;
    typedef Chrono::CpuTimePoint CpuTimePoint;
    typedef Chrono::CpuTimeDuration CpuTimeDuration;

    typedef std::forward_list<Key> TimerPath;

    Timer();
    Timer(const Key&);
    virtual ~Timer() = default;

    void tick(const Key& timerName);
    void tock(const Key& timerName);
    void tock();

    double lap(const Key& timerName) const ;
    double elapsed(const Key& timerName) const;

    std::ostream& report(const Key& timerName, std::ostream& os) const;
    std::ostream& report(std::ostream& os) const;
    std::ostream& report(const TimerPath&, std::ostream& os) const;

    void additionalTime(const TimerPath& path, WallTimeDuration additionalWall, CpuTimeDuration additionalCpu, int additionalTicks);
    TimerPath currentTimerNodePath() const;

    void reset();
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

std::ostream& operator<<(std::ostream& os, const Nextsim::Timer& tim);

#endif /* SRC_INCLUDE_TIMER_HPP */
