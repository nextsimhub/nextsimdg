/*
 * @file Chrono.hpp
 *
 * @date Oct 28, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef SRC_INCLUDE_CHRONO_HPP_
#define SRC_INCLUDE_CHRONO_HPP_

#include <chrono>
#include <ctime>

namespace Nextsim {

class Chrono {
public:

    typedef std::chrono::high_resolution_clock::time_point WallTimePoint;
    typedef std::chrono::high_resolution_clock::duration WallTimeDuration;

    typedef std::clock_t CpuTimePoint;
    typedef double CpuTimeDuration;

    Chrono()
        : m_wallTime(WallTimeDuration::zero())
        , m_cpuTime(0)
        , m_ticks(0)
        , m_running(false) {};
    ~Chrono() = default;

    inline WallTimePoint wallHack() const {return m_wallHack;};
    inline WallTimeDuration wallTime() const
    {
        return m_wallTime + (m_running ? wallTimeSinceHack() : WallTimeDuration::zero());
    };


    inline CpuTimePoint cpuHack() const {return m_cpuHack;};
    inline CpuTimeDuration cpuTime() const
    {
        return m_cpuTime + (m_running ? cpuTimeSinceHack() : 0);
    };

    inline int ticks() const {return m_ticks;};
    inline bool running() const {return m_running;};

    inline void start()
    {
        m_wallHack = std::chrono::high_resolution_clock::now();
        m_cpuHack = std::clock();
        ++m_ticks;
        m_running = true;
    };

    void stop()
    {
        m_wallTime += wallTimeSinceHack();
        m_cpuTime += cpuTimeSinceHack();
        m_running = false;
    };

    inline void extraCpuTime(const CpuTimeDuration& extraTime) {m_cpuTime += extraTime;};
    inline void extraWallTime(const WallTimeDuration& extraTime) {m_wallTime += extraTime;};
    inline void extraTicks(int extraTicks) {m_ticks += extraTicks;};

private:
    WallTimePoint m_wallHack;
    WallTimeDuration m_wallTime;

    CpuTimePoint m_cpuHack;
    CpuTimeDuration m_cpuTime;

    int m_ticks;

    bool m_running;

    CpuTimeDuration cpuTimeSinceHack() const
    {
        return static_cast<CpuTimeDuration>(std::clock() - m_cpuHack) / CLOCKS_PER_SEC;
    }

    WallTimeDuration wallTimeSinceHack() const
    {
        return std::chrono::duration_cast<WallTimeDuration>(
        std::chrono::high_resolution_clock::now() - m_wallHack);
    }

};

} /* namespace Nextsim */

#endif /* SRC_INCLUDE_CHRONO_HPP_ */
