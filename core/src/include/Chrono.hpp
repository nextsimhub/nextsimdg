/*!
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
/*!
 * @brief A class providing a timer.
 *
 * @details This class records execution wall and CPU time as well as a record
 * of the number of times that the chronometer has been started.
 */
class Chrono {
public:
    //! Type of a point in time for the wall clock.
    typedef std::chrono::high_resolution_clock::time_point WallTimePoint;
    //! Type of a time duration on the wall clock.
    typedef std::chrono::high_resolution_clock::duration WallTimeDuration;

    //! Type of a point in time for the CPU clock.
    typedef std::clock_t CpuTimePoint;
    //! Type of a time duration on the CPU clock.
    typedef double CpuTimeDuration;

    //! Default constructor of a chronometer, zeroing all of the counters.
    Chrono()
        : m_wallTime(WallTimeDuration::zero())
        , m_cpuTime(0)
        , m_ticks(0)
        , m_running(false) {};
    ~Chrono() = default;

    //! Returns the current time on the wall clock.
    inline WallTimePoint wallHack() const { return m_wallHack; };
    //! Returns the current cumulative wall clock time.
    inline WallTimeDuration wallTime() const
    {
        return m_wallTime + (m_running ? wallTimeSinceHack() : WallTimeDuration::zero());
    };

    //! Resets all of the chronometer counters.
    inline void reset()
    {
        m_wallTime = WallTimeDuration::zero();
        m_cpuTime = 0;
        m_ticks = 0;
        m_running = false;
    }

    //! Returns the current time on the CPU clock.
    inline CpuTimePoint cpuHack() const { return m_cpuHack; };
    //! Returns the current cumulative CPU clock timer.
    inline CpuTimeDuration cpuTime() const
    {
        return m_cpuTime + (m_running ? cpuTimeSinceHack() : 0);
    };

    //! Returns the current number of activation ticks.
    inline int ticks() const { return m_ticks; };
    //! Returns whether this chronometer is running.
    inline bool running() const { return m_running; };

    /*!
     * @brief Starts the timer.
     *
     * @details Starts the clock on both the wall and CPU clocks, increments
     * the number of activation ticks and sets the running flag.
     */
    inline void start()
    {
        m_wallHack = std::chrono::high_resolution_clock::now();
        m_cpuHack = std::clock();
        ++m_ticks;
        m_running = true;
    };

    /*!
     * @brief Stops the timer.
     *
     * @details Stops both the wall and CPU clocks, updates the cumulative
     * time for both clocks and unsets the running flag.
     */
    void stop()
    {
        m_wallTime += wallTimeSinceHack();
        m_cpuTime += cpuTimeSinceHack();
        m_running = false;
    };

    /*!
     * @brief Adds an externally determined increment to the CPU clock.
     *
     * @param extraTime the additional duration to be added to the CPU clock.
     */
    inline void extraCpuTime(const CpuTimeDuration& extraTime) { m_cpuTime += extraTime; };
    /*!
     * @brief Adds an externally determined increment to the wall clock.
     *
     * @param extraTime the additional duration to be added to the wall clock.
     */
    inline void extraWallTime(const WallTimeDuration& extraTime) { m_wallTime += extraTime; };
    /*!
     * @brief Adds an externally determined increment to the activation count.
     *
     * @param extraTicks the additional ticks to be added to the activation
     * count.
     */
    inline void extraTicks(int extraTicks) { m_ticks += extraTicks; };

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
