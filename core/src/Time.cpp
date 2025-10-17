/*!
 * @author  Tim Spain <timothy.spain@nersc.no>
 */

#include "include/Time.hpp"

#include <regex>
#include <sstream>
#include <stdexcept>

namespace Nextsim {

const std::string TimePoint::ymdFormat = "%Y-%m-%d";
const std::string TimePoint::doyFormat = "%Y-%j";
const std::string TimePoint::hmsFormat = "T%H:%M:%SZ";
const std::string TimePoint::ymdhmsFormat = ymdFormat + hmsFormat;
const std::string TimePoint::doyhmsFormat = doyFormat + hmsFormat;
std::array<bool, TimeOptions::COUNT> TimeOptions::m_opt = { false, false };

static const int minuteSeconds = 60;
static const int hourSeconds = minuteSeconds * 60;
static const int daySeconds = hourSeconds * 24;
static const int yearSeconds = daySeconds * 365;
static const int tmEpochYear = 1900;
static const int unixEpochYear = 1970;

Duration TimePoint::operator-(const TimePoint& a) const { return Duration(m_t - a.m_t); }
TimePoint& TimePoint::operator+=(const Duration& d)
{
    m_t += d.m_d;
    return *this;
}
TimePoint& TimePoint::operator-=(const Duration& d)
{
    m_t -= d.m_d;
    return *this;
}
TimePoint TimePoint::operator+(const Duration& d) const
{
    TimePoint t2(*this);
    return t2 += d;
}
bool TimePoint::operator<=(const TimePoint& a) const { return m_t <= a.m_t; }
bool TimePoint::operator<(const TimePoint& a) const { return m_t < a.m_t; }
bool TimePoint::operator>=(const TimePoint& a) const { return m_t >= a.m_t; }
bool TimePoint::operator>(const TimePoint& a) const { return m_t > a.m_t; }
bool TimePoint::operator==(const TimePoint& a) const { return m_t == a.m_t; }
bool TimePoint::operator!=(const TimePoint& a) const { return m_t != a.m_t; }

std::istream& TimePoint::parse(std::istream& is)
{
    auto fromTime = Clock::from_time_t(timeFromISO(is));
    m_t = fromTime;
    return is;
}

TimePoint& TimePoint::parse(const std::string& str)
{
    std::stringstream is(str);
    parse(is);
    return *this;
}

std::ostream& TimePoint::format(std::ostream& os, std::string formatStr) const
{
    // Temporary conversion from int to system_clock
    auto tt = Clock::to_time_t(m_t);
    os << std::put_time(std::gmtime(&tt), formatStr.c_str());
    return os;
}

std::string TimePoint::format(std::string formatStr) const
{
    std::stringstream ss;
    format(ss, formatStr);
    return ss.str();
}

std::tm& tmDoy(std::tm& tm)
{
    int common0th[] = { -1, 30, 58, 89, 119, 150, 180, 211, 242, 272, 303, 333 };
    int leap0th[] =   { -1, 30, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334 };
    int trueYear = tm.tm_year + tmEpochYear;
    bool isLeap = ((trueYear % 4 == 0) && (trueYear % 100 != 0)) || (trueYear % 400 == 0);
    int* zerothArray = isLeap ? leap0th : common0th;
    tm.tm_yday = zerothArray[tm.tm_mon] + tm.tm_mday;
    return tm;
}

const int TimePoint::year() const
{
    return gmtime()->tm_year + tmEpochYear;
}
const size_t TimePoint::month() const
{
    return gmtime()->tm_mon + 1;
}
const size_t TimePoint::day() const
{
    return gmtime()->tm_mday;
}
const size_t TimePoint::doy() const
{
    return gmtime()->tm_yday + 1;
}
const size_t TimePoint::hour() const
{
    return gmtime()->tm_hour;
}
const size_t TimePoint::minute() const
{
    return gmtime()->tm_min;
}
const double TimePoint::second() const
{
    return gmtime()->tm_sec;
}

std::time_t mkgmtime(std::tm* tm, bool recalculateDoy)
{
    if (recalculateDoy)
        tmDoy(*tm);
    std::time_t sum = tm->tm_sec;
    sum += tm->tm_min * minuteSeconds;
    sum += tm->tm_hour * hourSeconds;
    sum += tm->tm_yday * daySeconds;
    const int year = tmEpochYear + tm->tm_year;

    static const size_t yOffset = 399;
    static const size_t oEpochYear = unixEpochYear + yOffset;
    static const size_t epochDays = unixEpochYear * 365 + oEpochYear/4 - oEpochYear/100 + oEpochYear/400;
    const size_t oYear = year + yOffset; // Offset the year to make the leap years tick over correctyl
    const size_t startYearDays = year * 365 + oYear / 4 - oYear / 100 + oYear / 400;
    sum += (startYearDays  - epochDays) * daySeconds;

    return sum;
}

bool isYMDFormat(const std::string& iso)
{
    const std::regex ymd("^\\d+-\\d+-\\d+($|T)"); // Search for the month
    return std::regex_search(iso, ymd);
}

bool isDOYFormat(const std::string& iso)
{
    const std::regex doy("^\\d+-\\d+($|T)"); // Search for the day of year

    bool isDOY = std::regex_search(iso, doy);

    if (!isYMDFormat(iso) && !isDOY)
        throw std::invalid_argument("Unrecognized date format: " + iso);

    if (TimeOptions::useDOY() && !isDOY)
        throw std::invalid_argument("Inconsistent date format: " + iso
            + " with useDOY = " + (TimeOptions::useDOY() ? "true" : "false"));

    return isDOY;
}

std::vector<std::string> splitString(const std::string& str, char delim)
{
    std::stringstream ss(str);
    std::string s;
    std::vector<std::string> vs;
    while (std::getline(ss, s, delim)) {
        vs.push_back(s);
    }
    return vs;
}

std::string addTimeLeadingZeros(const std::string& in)
{
    // Linux stdlibc++ doesn't like time components without leading zeros
    std::vector<std::string> dateTime = splitString(in, 'T');
    std::string out = dateTime[0];
    if (dateTime.size() > 1) {
        std::vector<std::string> hms = splitString(dateTime[1], ':');
        std::stringstream timeStream;
        timeStream << std::setfill('0');
        for (size_t i = 0; i < hms.size(); ++i) {
            timeStream << std::setw(2) << std::stoi(hms[i]);
            if (i != hms.size() - 1)
                timeStream << ':';
        }
        out += "T" + timeStream.str();
    }

    return out;
}
std::tm getTimTime(const std::string& in, bool isDOY)
{
    std::tm tm;
    // Reset the time values
    tm.tm_hour = 0;
    tm.tm_min = 0;
    tm.tm_sec = 0;

    std::string iso = addTimeLeadingZeros(in);
    if (isDOY) {
        // Parse the string manually to deal with broken %j format in libstdc++
        // Split the string into date and time and prepare to parse the time portion later
        std::vector<std::string> dateTime = splitString(iso, 'T');
        if (dateTime.size() > 1) {
            std::stringstream timeStream("T" + dateTime[1]);
            timeStream >> std::get_time(&tm, TimePoint::hmsFormat.c_str());
        }
        // Parse the date portion by splitting on the (lone) hyphen
        std::vector<std::string> yearDoy = splitString(dateTime[0], '-');
        tm.tm_year = std::stoi(yearDoy[0]) - tmEpochYear;
        tm.tm_yday = std::stoi(yearDoy[1]) - 1;
    } else {
        std::stringstream(iso) >> std::get_time(&tm, TimePoint::ymdhmsFormat.c_str());
    }
    return tm;
}

std::time_t timeFromISO(const std::string& iso)
{
    bool isDOY = isDOYFormat(iso);
    std::tm tm = getTimTime(iso, isDOY);
    return mkgmtime(&tm, !isDOY);
}

std::time_t timeFromISO(std::istream& is)
{
    std::string iso;
    is >> iso;
    return timeFromISO(iso);
}

Duration durationFromISO(const std::string& iso, int sign = +1)
{
    if (isYMDFormat(iso)) {
        throw std::invalid_argument(
            "Duration does not accept months as they have arbitrary length");
    }
    bool isDOY = isDOYFormat(iso);
    std::tm tm = getTimTime(iso, isDOY);
    // Make up the time duration, analogously to mkgmtime()
    size_t sum = tm.tm_sec;
    sum += tm.tm_min * minuteSeconds;
    sum += tm.tm_hour * hourSeconds;
    if (isDOY) {
        sum += (tm.tm_yday + 1) * daySeconds;
    } else {
        // 30 day months until real calendar intervals are implemented
        sum += tm.tm_mon * 30 * daySeconds;
        sum += tm.tm_mday * daySeconds;
    }
    sum += (tmEpochYear + tm.tm_year) * yearSeconds;
    Duration::Basis dura(std::chrono::seconds(sign * sum));
    return Duration(dura);
}

Duration durationFromISO(std::istream& is, int sign = +1)
{
    std::string iso;
    is >> iso;
    return durationFromISO(iso, sign);
}

std::istream& Duration::parse(std::istream& is)
{
    // read the first character, check it is the ISO standard P
    char possibleP;
    is >> possibleP;
    if (possibleP != 'P') {
        std::string restOf;
        is >> restOf;
        restOf = possibleP + restOf;
        // If the remaining string is a valid double, then interpret that as a
        // duration in seconds, else throw an invalid argument exception.

        // Double regex courtesy of https://stackoverflow.com/a/56502134
        std::regex rx(R"(^([+-]?(?:[[:d:]]+\.?|[[:d:]]*\.[[:d:]]+))(?:[Ee][+-]?[[:d:]]+)?$)");
        bool isYMD = std::regex_search(restOf, rx);
        if (!isYMD)
            throw std::invalid_argument(
                "The duration should be an ISO 8601 duration (P…) or a number of seconds. Got: "
                + restOf);
        double sec = std::stod(restOf);
        // Assign the seconds value to the Duration
        setDurationSeconds(sec);
        return is;
    }

    // Peek at the next character, to see if it is a -
    bool isNegative = (is.peek() == '-');
    if (isNegative) {
        // pop the negative sign, then parse the rest
        char sign;
        is >> sign;
    }
    *this = durationFromISO(is, isNegative ? -1 : 1);
    // Temporary conversion from system_clock to int
    return is;
}

Duration::Duration(const std::string& str) { this->parse(str); }

Duration::Duration(double seconds) { setDurationSeconds(seconds); }

void Duration::setDurationSeconds(double secs)
{
    std::chrono::duration<double> sec(secs);
    m_d = std::chrono::duration_cast<Basis>(sec);
}

TimePoint Duration::operator+(const TimePoint& t) const { return t + *this; }

std::tm* TimePoint::gmtime() const
{
    auto tt = Clock::to_time_t(m_t);
    return std::gmtime(&tt);
}
}
