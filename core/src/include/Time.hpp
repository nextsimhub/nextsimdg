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

class TimePointImpl {
public:
    TimePointImpl();
    TimePointImpl(const std::string&);
    TimePointImpl(const TimePoint&);
    TimePointImpl(const TimePoint&, const Duration&);

    Duration operator-(const TimePoint& a) { return m_t - a; }
    TimePoint& operator+=(const Duration& d)
    {
        m_t += d;
        return *this;
    }
    TimePoint operator+(const Duration& d) { return m_t + d; }

    bool operator<=(const TimePoint& a) { return m_t <= a; }
    bool operator<(const TimePoint& a) { return m_t < a; }
    bool operator>=(const TimePoint& a) { return m_t >= a; }
    bool operator>(const TimePoint& a) { return m_t > a; }
    bool operator==(const TimePoint& a) { return m_t == a; }
    bool operator!=(const TimePoint& a) { return m_t != a; }
private:
    TimePoint m_t; //FIXME: Once implemented, change this to a chrono::time_point and the typedef to point here.
};

class DurationImpl {
public:
    DurationImpl();
    DurationImpl(const std::string&);
    DurationImpl(const Duration&);

    TimePoint operator+(const TimePoint& t) const { return t + *this; }

    Duration& operator+=(const Duration& a)
    {
        m_d += a;
        return *this;
    }
    Duration& operator-=(const Duration& a)
    {
        m_d -= a;
        return *this;
    }

    Duration& operator*=(double a)
    {
        m_d *= a;
        return *this;
    }
    Duration& operator/=(double a)
    {
        m_d /= a;
        return *this;
    }

    Duration operator+(const DurationImpl& a) const
    {
        Duration d = m_d;
        return d + a.m_d;
    }
    Duration operator-(const DurationImpl& a) const
    {
        Duration d = m_d;
        return d - a.m_d;
    }

    double seconds() const { return m_d; }
private:
    Duration m_d; //FIXME: Once implemented, change this to a chrono::duration and the typedef to point here.
};

inline double operator*(double a, const DurationImpl& b) { return a * b.seconds(); }
inline double operator/(double a, const DurationImpl& b) { return a / b.seconds(); }
inline double operator*(const DurationImpl& a, double b) { return b * a; }
inline double operator/(const DurationImpl& a, double b) { return a.seconds() / b; }
}; // namespace Nextsim

#endif /* TIME_HPP */
