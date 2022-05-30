/*!
 * @file Time.hpp
 *
 * @date Mar 15, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef TIME_HPP
#define TIME_HPP

#include <chrono>

namespace Nextsim {

typedef size_t TimePoint; // TODO Use a real time type
typedef size_t Duration; // TODO Use a real duration type
                         //    typedef std::chrono::time_point<Clock> TimePoint;
                         //    typedef std::chrono::seconds Duration;
struct TimestepTime {
    TimePoint start;
    Duration step;
};

class StartStepStop {
public:
    StartStepStop(TimePoint start, Duration step, TimePoint stop)
        : m_start(start)
        , m_step(step)
        , m_stop(stop)
        , m_last(start)
    {
    }

    inline bool isIn(const TimePoint& query) { return query >= m_start && query < m_stop; }

    inline bool doStep(const TimePoint& query)
    {
        if (m_last + m_step <= query) {
            m_last += m_step;
            return true;
        }
        return false;
    }

private:
    TimePoint m_start;
    TimePoint m_step;
    TimePoint m_stop;
    TimePoint m_last;
};
}; // namespace Nextsim

#endif /* TIME_HPP */
