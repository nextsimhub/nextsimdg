/*!
 * @file Time.hpp
 *
 * @date Mar 15, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef TIME_HPP
#define TIME_HPP

#include <array>
#include <chrono>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <sstream>

namespace Nextsim {

class TimeOptions {
    enum TimeOption { USE_DOY, DAYS_360, COUNT };

public:
    static bool& useDOY() { return m_opt[USE_DOY]; }

private:
    static std::array<bool, COUNT> m_opt;
};

typedef std::chrono::system_clock SystemClock;
typedef SystemClock::duration SystemDuration;

class TimePoint;

class Duration {
public:
    typedef SystemDuration Basis;
    Duration()
        : m_d()
    {
    }
    Duration(const std::string& str) { this->parse(str); }

    TimePoint operator+(const TimePoint& t) const;

    Duration& operator+=(const Duration& a)
    {
        m_d += a.m_d;
        return *this;
    }
    Duration& operator-=(const Duration& a)
    {
        m_d -= a.m_d;
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

    Duration operator+(const Duration& a) const { return Duration(m_d + a.m_d); }
    Duration operator-(const Duration& a) const { return Duration(m_d - a.m_d); }

    double seconds() const { return std::chrono::duration_cast<std::chrono::seconds>(m_d).count(); }

    std::istream& parse(std::istream& is);

    Duration& parse(const std::string& str)
    {
        std::stringstream is(str);
        parse(is);
        return *this;
    }

    std::ostream& format(std::ostream& os) const { return os << seconds(); }

    std::string format() const
    {
        std::stringstream ss;
        format(ss);
        return ss.str();
    }

    friend TimePoint;
    friend Duration durationFromISO(const std::string&, int);

private:
    Duration(const Basis& d)
        : m_d(d)
    {
    }
    Basis m_d;
};

/*!
 * @brief converts the std::tm struct to a std::time_t value assuming UTC.
 *
 * @details A UTC-fixed version of the standard library function mktime().
 *          Calculates the duration between the Unix epoch and the given
 *          std::tm time and directly assigns that value.
 *
 * @param time Pointer to the tm structure to interpret.
 * @param recalculateDoy Flag argument for whether to calculate the day of the
 *          year from the year month and day values.
 */
std::time_t mkgmtime(std::tm* time, bool recalculateDoy = true);

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

class TimePoint {
public:
    typedef SystemClock Clock;
    typedef std::chrono::time_point<Clock, Duration::Basis> Basis;

    TimePoint()
        : m_t() {};
    TimePoint(const std::string& str) { this->parse(str); }
    TimePoint(const TimePoint&, const Duration&);

    Duration operator-(const TimePoint& a) const { return Duration(m_t - a.m_t); }
    TimePoint& operator+=(const Duration& d)
    {
        m_t += d.m_d;
        return *this;
    }
    TimePoint& operator-=(const Duration& d)
    {
        m_t -= d.m_d;
        return *this;
    }
    TimePoint operator+(const Duration& d) const
    {
        TimePoint t2(*this);
        return t2 += d;
    }

    bool operator<=(const TimePoint& a) const { return m_t <= a.m_t; }
    bool operator<(const TimePoint& a) const { return m_t < a.m_t; }
    bool operator>=(const TimePoint& a) const { return m_t >= a.m_t; }
    bool operator>(const TimePoint& a) const { return m_t > a.m_t; }
    bool operator==(const TimePoint& a) const { return m_t == a.m_t; }
    bool operator!=(const TimePoint& a) const { return m_t != a.m_t; }

    std::istream& parse(std::istream& is)
    {
        auto fromTime = Clock::from_time_t(timeFromISO(is));
        m_t = fromTime;
        return is;
    }

    TimePoint& parse(const std::string& str)
    {
        std::stringstream is(str);
        parse(is);
        return *this;
    }

    std::ostream& format(std::ostream& os) const
    {
        // Temporary conversion from int to system_clock
        auto tt = Clock::to_time_t(m_t);
        os << std::put_time(std::gmtime(&tt), ymdhmsFormat.c_str());
        return os;
    }

    std::string format() const
    {
        std::stringstream ss;
        format(ss);
        return ss.str();
    }
    // FIXME Remove me
    Basis& getTime() { return m_t; }

    std::tm* gmtime() const;

    static const std::string ymdFormat;
    static const std::string doyFormat;
    static const std::string ymdhmsFormat;
    static const std::string doyhmsFormat;
    static const std::string hmsFormat;

    friend Duration;

private:
    TimePoint(const Basis& t)
        : m_t(t)
    {
    }
    Basis m_t;
};

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
    Duration m_step;
    TimePoint m_stop;
    TimePoint m_last;
};

inline double operator*(double a, const Duration& b) { return a * b.seconds(); }
inline double operator/(double a, const Duration& b) { return a / b.seconds(); }
inline double operator*(const Duration& a, double b) { return b * a; }
inline double operator/(const Duration& a, double b) { return a.seconds() / b; }

inline std::istream& operator>>(std::istream& is, TimePoint& tp) { return tp.parse(is); }
inline std::istream& operator>>(std::istream& is, Duration& dur) { return dur.parse(is); }
inline std::ostream& operator<<(std::ostream& os, const TimePoint& tp) { return tp.format(os); }
inline std::ostream& operator<<(std::ostream& os, const Duration& dur) { return dur.format(os); }

}; // namespace Nextsim

#endif /* TIME_HPP */
