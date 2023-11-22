/*!
 * @file Time.cpp
 *
 * @date Jun 7, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
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

std::tm& tmDoy(std::tm& tm)
{
    int month0th[] = { 0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334 };
    bool isLeap = ((tm.tm_year % 4 == 0) && (tm.tm_year % 100 != 0)) || (tm.tm_year % 400 == 0);
    int bissextile = (isLeap && tm.tm_mon >= 2) ? 1 : 0;
    tm.tm_yday = month0th[tm.tm_mon] + tm.tm_mday + bissextile;
    return tm;
}

std::time_t mkgmtime(std::tm* tm, bool recalculateDoy)
{
    if (recalculateDoy)
        tmDoy(*tm);
    std::time_t sum = tm->tm_sec;
    sum += tm->tm_min * minuteSeconds;
    sum += tm->tm_hour * hourSeconds;
    sum += (tm->tm_yday - 1) * daySeconds;
    int year = tmEpochYear + tm->tm_year;
    std::time_t unixYear = year - unixEpochYear;
    sum += unixYear * yearSeconds;

    // Handle the effect of leap days on the first day of the year. Proleptic Gregorian.
    int julianLeapDays = (year - 1) / 4 - (unixEpochYear - 1) / 4;
    sum += julianLeapDays * daySeconds;
    // Skipped Gregorian leap days
    sum += (julianGregorianShiftDays(year) - julianGregorianShiftDays(unixEpochYear)) * daySeconds;

    return sum;
}

int julianGregorianShiftDays(int year)
{
    // Note the subtraction of 1 year, as xx00 behaves like (xx-1)99 and not
    // necessarily like xx01
    int centurySinceGreg = (year - 1) / 100 - 15;
    int leaps = (3 * centurySinceGreg) / 4 + 10;
    return -leaps;
}

bool isDOYFormat(const std::string& iso)
{
    const std::regex ymd("^\\d+-\\d+-\\d+($|T)"); // Search for the month
    const std::regex doy("^\\d+-\\d+($|T)"); // Search for the day of year

    bool isYMD = std::regex_search(iso, ymd);
    bool isDOY = std::regex_search(iso, doy);

    if (!isYMD && !isDOY)
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
        tm.tm_yday = std::stoi(yearDoy[1]);
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
    bool isDOY = isDOYFormat(iso);
    std::tm tm = getTimTime(iso, isDOY);
    // Make up the time duration, analogously to mkgmtime()
    size_t sum = tm.tm_sec;
    sum += tm.tm_min * minuteSeconds;
    sum += tm.tm_hour * hourSeconds;
    if (isDOY) {
        sum += tm.tm_yday * daySeconds;
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
