/*!
 * @file Time.cpp
 *
 * @date Jun 7, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/Time.hpp"

#include <regex>
#include <sstream>

namespace Nextsim {

const std::string TimePointImpl::ymdFormat = "%Y-%m-%d";
const std::string TimePointImpl::doyFormat = "%Y-%j";
const std::string TimePointImpl::hmsFormat = "T%H:%M:%SZ";
const std::string TimePointImpl::ymdhmsFormat = ymdFormat + hmsFormat;
const std::string TimePointImpl::doyhmsFormat = doyFormat + hmsFormat;
std::array<bool, TimeOptions::COUNT> TimeOptions::m_opt = { false, false };

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
    const int minuteSeconds = 60;
    const int hourSeconds = minuteSeconds * 60;
    const int daySeconds = hourSeconds * 24;
    const int yearSeconds = daySeconds * 365;
    const int tmEpochYear = 1900;
    const int unixEpochYear = 1970;

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

std::time_t timeFromISO(const std::string& iso)
{
    const std::regex ymd("-\\d\\d-"); // Search for the month (exactly two digits)
    const std::regex doy("-\\d\\d\\d"); // Search for the day of year (exactly three digits)

    bool isYMD = std::regex_search(iso, ymd);
    bool isDOY = std::regex_search(iso, doy);

    if (!isYMD && !isDOY)
        throw std::invalid_argument("Unrecognized date format: " + iso);

    if (TimeOptions::useDOY() && !isDOY)
        throw std::invalid_argument("Inconsistent date format: " + iso
                + " with useDOY = " + (TimeOptions::useDOY() ? "true" : "false"));

    std::tm tm;
    // Reset the time values
    tm.tm_hour = 0;
    tm.tm_min = 0;
    tm.tm_sec = 0;

    const char* formatCStr
        = (isYMD) ? TimePointImpl::ymdhmsFormat.c_str() : TimePointImpl::doyhmsFormat.c_str();

    std::stringstream isoStream(iso);
    isoStream >> std::get_time(&tm, formatCStr);
    return mkgmtime(&tm, isYMD);
}

std::time_t timeFromISO(std::istream& is)
{
    std::string iso;
    is >> iso;
    return timeFromISO(iso);
}

}
