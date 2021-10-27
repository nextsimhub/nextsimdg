/*!
 * @file Timer.hpp
 * @date 11 Aug 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef SRC_INCLUDE_TIMER_HPP
#define SRC_INCLUDE_TIMER_HPP

#include "TimerBase.hpp"

#include <chrono>
#include <ctime>
#include <forward_list>
#include <map>
#include <ostream>
#include <set>
#include <string>

namespace Nextsim {

class Timer : public TimerBase {
public:
    typedef std::chrono::high_resolution_clock::time_point WallTimePoint;
    typedef std::chrono::high_resolution_clock::duration WallTimeDuration;

    typedef std::clock_t CpuTimePoint;
    typedef double CpuTimeDuration;

    typedef std::forward_list<Key> TimerPath;

    Timer();
    Timer(const Key&);
    virtual ~Timer() = default;

    void tick(const Key& timerName) override;
    void tock(const Key& timerName) override;

    double lap(const Key& timerName) const override;
    double elapsed(const Key& timerName) const override;

    std::ostream& report(const Key& timerName, std::ostream& os) const override;
    std::ostream& report(std::ostream& os) const override;
    std::ostream& report(const TimerPath&, std::ostream& os) const;

    void additionalTime(const TimerPath& path, WallTimeDuration additionalWall, CpuTimeDuration additionalCpu, int additionalTicks);
    TimerPath currentTimerNodePath() const;

    static Timer main;

private:

    struct TimerNode {
        TimerNode();
        Key name;

        WallTimePoint wallHack;
        WallTimeDuration wallTime;

        CpuTimePoint cpuHack;
        CpuTimeDuration cpuTime;

        int ticks;

        bool running;

        std::map<Key, TimerNode> childNodes;
        TimerNode* parent;

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
