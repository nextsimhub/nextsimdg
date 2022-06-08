/*!
 * @file Time.hpp
 *
 * @date Mar 15, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef TIME_HPP
#define TIME_HPP

#include <chrono>
#include <ctime>
#include <iomanip>
#include <iostream>

namespace Nextsim {

typedef size_t TimePoint; // TODO Use a real time type
typedef int Duration; // TODO Use a real duration type
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

/*!
 * @brief converts the std::tm struct to a std::time_t value assuming UTC.
 *
 * @details A UTC-fixed version of the standard library function mktime().
 *          Calculates the duration between the Unix epoch and the given
 *          std::tm time and directly assigns that value.
 *
 * @param time Pointer to the tm structure to interpret.
 */
std::time_t mkgmtime(std::tm* time);

/*!
 * @brief Returns the number of days between the Julian and Gregorian years.
 *
 * @details Returns the number of days between 1st January of the given year in
 * the Julian and Gregorian calendars.
 */
int julianGregorianShiftDays(int year);

/*!
 * @brief Returns a std::time_t from a given ISO date string of a specific format.
 *
 * @param iso A std::string in the ISO YYYY-MM-DDThh:mm:ss format.
 */
std::time_t timeFromISO(const std::string& iso);
/*!
 * @brief Returns a std::time_t from a given ISO date string of a specific format.
 *
 * @param iso A std::istream containing an ISO YYYY-MM-DDThh:mm:ss format date.
 */
std::time_t timeFromISO(std::istream& is);
class TimePointImpl {
public:
    TimePointImpl()
        : m_t() {};
    TimePointImpl(const std::string&);
    TimePointImpl(const TimePoint&);
    TimePointImpl(const TimePoint&, const Duration&);

    Duration operator-(const TimePoint& a) { return m_t - a; }
    TimePointImpl& operator+=(const Duration& d)
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

    std::istream& parse(std::istream& is)
    {
        auto interTime = std::chrono::system_clock::from_time_t(timeFromISO(is));//mkgmtime(&tm));
        // Temporary conversion from system_clock to int
        m_t = interTime.time_since_epoch().count() / 1000000;
        return is;
    }

    std::ostream& format(std::ostream& os)
    {
        // Temporary conversion from int to system_clock
        std::chrono::duration<int> sinceEpoch(m_t);
        std::chrono::time_point<std::chrono::system_clock, std::chrono::duration<int>> sysTime(
            sinceEpoch);
        auto tTime = std::chrono::system_clock::to_time_t(sysTime);
        os << std::put_time(std::gmtime(&tTime), isoFormatString);
        return os;
    }

    // FIXME Remove me
    TimePoint getTime() { return m_t; }

    static const char isoFormatString[];
private:
    TimePoint m_t; // FIXME: Once implemented, change this to a chrono::time_point and the typedef
                   // to point here.
};

class DurationImpl {
public:
    DurationImpl();
    DurationImpl(const std::string&);
    DurationImpl(const Duration&);

    TimePoint operator+(const TimePoint& t) const { return t + this->m_d; }

    DurationImpl& operator+=(const Duration& a)
    {
        m_d += a;
        return *this;
    }
    DurationImpl& operator-=(const Duration& a)
    {
        m_d -= a;
        return *this;
    }

    DurationImpl& operator*=(double a)
    {
        m_d *= a;
        return *this;
    }
    DurationImpl& operator/=(double a)
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
    Duration m_d; // FIXME: Once implemented, change this to a chrono::duration and the typedef to
                  // point here.
};

inline double operator*(double a, const DurationImpl& b) { return a * b.seconds(); }
inline double operator/(double a, const DurationImpl& b) { return a / b.seconds(); }
inline double operator*(const DurationImpl& a, double b) { return b * a; }
inline double operator/(const DurationImpl& a, double b) { return a.seconds() / b; }

}; // namespace Nextsim

#endif /* TIME_HPP */
