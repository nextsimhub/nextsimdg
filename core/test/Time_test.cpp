/*!
 * @author  Tim Spain <timothy.spain@nersc.no>
 */

#include "include/Time.hpp"

#include <sstream>

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

namespace Nextsim {

void testTimePointTime(int year, size_t month, size_t day, size_t hour, size_t minute, double second)
{
    std::string timeStr = std::to_string(year)+
            "-"+std::to_string(month)+
            "-"+std::to_string(day)+
            "T"+std::to_string(hour)+
            ":"+std::to_string(minute)+
            ":"+std::to_string(second);
    TimePoint testTime(timeStr);
    REQUIRE_MESSAGE(testTime.year() == year, "for TimePoint representing ", timeStr);
    REQUIRE_MESSAGE(testTime.month() == month, "for TimePoint representing ", timeStr);
    REQUIRE_MESSAGE(testTime.day() == day, "for TimePoint representing ", timeStr);
    REQUIRE_MESSAGE(testTime.hour() == hour, "for TimePoint representing ", timeStr);
    REQUIRE_MESSAGE(testTime.minute() == minute, "for TimePoint representing ", timeStr);
    REQUIRE_MESSAGE(testTime.second() == second, "for TimePoint representing ", timeStr);
}

TEST_SUITE_BEGIN("Time");
TEST_CASE("TimePoint parsing and formating")
{
    // Time with explicit timezone marker so the initial and final strings match.
    std::stringstream is("2022-06-07T14:16:01Z");

    TimePoint tp;
    tp.parse(is);

    std::stringstream os;
    tp.format(os);

    REQUIRE(is.str() == os.str());

    // Directly with strings
    TimePoint tq;
    std::string rightNow = "2022-06-10T09:33:21Z";
    tq.parse(rightNow);
    REQUIRE(tq.format() == rightNow);
}

TEST_CASE("tmDoy")
{
    const char* iso = TimePoint::ymdhmsFormat.c_str();
    std::tm& tmDoy(std::tm& tm);
    std::tm tm;

    // 1995: Common year, yday(28th Feb) = 58, yday(1st Mar) = 59
    std::stringstream("1995-02-28T00:00:00") >> std::get_time(&tm, iso);
    tm.tm_yday = 0; // Reset the day of the year
    tmDoy(tm);
    REQUIRE(tm.tm_yday == 58);

    std::stringstream("1995-03-01T00:00:00") >> std::get_time(&tm, iso);
    tm.tm_yday = 0; // Reset the day of the year
    tmDoy(tm);
    REQUIRE(tm.tm_yday == 59);

    // 1996: Leap year, yday(28th Feb) = 58, yday(29th Feb) = 59, yday(1st Mar) = 60
    std::stringstream("1996-02-28T00:00:00") >> std::get_time(&tm, iso);
    tm.tm_yday = 0; // Reset the day of the year
    tmDoy(tm);
    REQUIRE(tm.tm_yday == 58);

    std::stringstream("1996-02-29T00:00:00") >> std::get_time(&tm, iso);
    tm.tm_yday = 0; // Reset the day of the year
    tmDoy(tm);
    REQUIRE(tm.tm_yday == 59);

    std::stringstream("1996-03-01T00:00:00") >> std::get_time(&tm, iso);
    tm.tm_yday = 0; // Reset the day of the year
    tmDoy(tm);
    REQUIRE(tm.tm_yday == 60);

    // 2000: Leap year, yday(28th Feb) = 58, yday(29th Feb) = 59, yday(1st Mar) = 60
    std::stringstream("2000-02-28T00:00:00") >> std::get_time(&tm, iso);
    tm.tm_yday = 0; // Reset the day of the year
    tmDoy(tm);
    REQUIRE(tm.tm_yday == 58);

    std::stringstream("2000-02-29T00:00:00") >> std::get_time(&tm, iso);
    tm.tm_yday = 0; // Reset the day of the year
    tmDoy(tm);
    REQUIRE(tm.tm_yday == 59);

    std::stringstream("2000-03-01T00:00:00") >> std::get_time(&tm, iso);
    tm.tm_yday = 0; // Reset the day of the year
    tmDoy(tm);
    REQUIRE(tm.tm_yday == 60);

    // 1900: Common year, yday(28th Feb) = 58, yday(1st Mar) = 59
    std::stringstream("1900-02-28T00:00:00") >> std::get_time(&tm, iso);
    tm.tm_yday = 0; // Reset the day of the year
    tmDoy(tm);
    REQUIRE(tm.tm_yday == 58);

    std::stringstream("1900-03-01T00:00:00") >> std::get_time(&tm, iso);
    tm.tm_yday = 0; // Reset the day of the year
    tmDoy(tm);
    REQUIRE(tm.tm_yday == 59);
}

TEST_CASE("mkgmtime")
{
    const int days = 24 * 60 * 60;

    const char* iso = TimePoint::ymdhmsFormat.c_str();

    // Leap days
    std::tm sextile;
    std::tm calends;
    // 1999: Ordinary year
    std::stringstream("1999-02-28T00:00:00") >> std::get_time(&sextile, iso);
    std::stringstream("1999-03-01T00:00:00") >> std::get_time(&calends, iso);
    double dayDiff = std::difftime(mkgmtime(&calends), mkgmtime(&sextile));
    REQUIRE(dayDiff == days);

    // 2000: Ordinary year
    std::stringstream("2000-02-28T00:00:00") >> std::get_time(&sextile, iso);
    std::stringstream("2000-03-01T00:00:00") >> std::get_time(&calends, iso);
    dayDiff = std::difftime(mkgmtime(&calends), mkgmtime(&sextile));
    REQUIRE(dayDiff == 2*days);

    std::stringstream ss("1970-01-01T00:00:00");
    std::tm epoch_tm;
    ss >> std::get_time(&epoch_tm, iso);

    std::time_t epoch_time = mkgmtime(&epoch_tm);

    REQUIRE(epoch_time == 0);

    std::stringstream ss1("1971-01-01T00:00:00");
    std::tm tm;
    //    tmZero(tm);
    ss1 >> std::get_time(&tm, iso);
    REQUIRE(mkgmtime(&tm) == 365 * days);

    std::tm tow_tm;
    std::stringstream stow("2022-06-07T15:37:45"); // UTC
    stow >> std::get_time(&tow_tm, iso);
    std::time_t timeOfWriting = mkgmtime(&tow_tm);

    REQUIRE(timeOfWriting == 1654616265);

    // Leap Day William
    std::stringstream("1980-02-28T00:00:00") >> std::get_time(&tm, iso);
    std::time_t before = mkgmtime(&tm);

    std::stringstream("1980-03-01T00:00:00") >> std::get_time(&tm, iso);
    std::time_t after = mkgmtime(&tm);

    REQUIRE((after - before) == 2 * days);
    std::stringstream("1980-02-13T00:00:00") >> std::get_time(&tm, iso);
    REQUIRE(mkgmtime(&tm) == 319248000);

    // Gregorian non-leap day
    std::stringstream("2100-02-28T00:00:00") >> std::get_time(&tm, iso);
    before = mkgmtime(&tm);

    std::stringstream("2100-03-01T00:00:00") >> std::get_time(&tm, iso);
    after = mkgmtime(&tm);

    REQUIRE((after - before) == 1 * days);
    //    REQUIRE(std::mktime(&tm) > 0); // How does linux return -ve here?
    REQUIRE(before == 4107456000);
    REQUIRE(after == 4107542400);

    // Turn of the 22nd century
    std::stringstream("2099-12-31T00:00:00") >> std::get_time(&tm, iso);
    REQUIRE(mkgmtime(&tm) == 4102358400);
    std::stringstream("2100-01-01T00:00:00") >> std::get_time(&tm, iso);
    REQUIRE(mkgmtime(&tm) == 4102444800);
    std::stringstream("2101-01-01T00:00:00") >> std::get_time(&tm, iso);
    REQUIRE(mkgmtime(&tm) == 4133980800);
    std::stringstream("2101-01-01T00:00:00") >> std::get_time(&tm, iso);
    after = mkgmtime(&tm);
    std::stringstream("2100-01-01T00:00:00") >> std::get_time(&tm, iso);
    before = mkgmtime(&tm);
    REQUIRE(after - before == 365 * days);
}

TEST_CASE("timeFromISO")
{
    const int days = 24 * 60 * 60;

    // Not a date, should throw
    REQUIRE_THROWS(timeFromISO("fhqwhgads"));

    // Default format
    REQUIRE(timeFromISO("1970-01-03T12:00:00Z") == 2 * days + (days / 2));
    // No Zulu time marker
    REQUIRE(timeFromISO("1970-02-01T06:00:00") == 31 * days + (days / 4));
    // No time value (defaults to 00UT)
    REQUIRE(timeFromISO("1971-01-02") == 366 * days);

    // YYYY-DDD parsing
    REQUIRE(timeFromISO("1971-003T00:00:00") == 367 * days);
    REQUIRE(timeFromISO("1971-031") == 395 * days);

    // Set DoY only parsing, try reading YMD
    bool previousDOY = TimeOptions::useDOY();
    TimeOptions::useDOY() = true;
    REQUIRE_THROWS(timeFromISO("1971-01-02"));
    TimeOptions::useDOY() = previousDOY;
}

TEST_CASE("TimePoints")
{
    TimePoint tp;
    TimePoint tq;

    // Comparison operators
    REQUIRE(tp.parse("1980-07-30") == tq.parse("1980-212"));
    REQUIRE(tp.parse("1980-07-30") > tq.parse("1980-07-29"));
    REQUIRE(tp.parse("1980-07-30") >= tq.parse("1980-212"));
    REQUIRE(tp.parse("1980-07-30") >= tq.parse("1980-07-29"));
    REQUIRE(tp.parse("1980-07-30") < tq.parse("1980-07-31"));
    REQUIRE(tp.parse("1980-07-30") <= tq.parse("1980-212"));
    REQUIRE(tp.parse("1980-07-30") <= tq.parse("1980-07-31"));
    REQUIRE(tp.parse("1980-07-30") != tq.parse("1980-07-29"));
}

TEST_CASE("Durations")
{
    const int days = 24 * 60 * 60;

    Duration dur;
    // Basic values
    REQUIRE_THROWS(dur.parse("0-0-0T0:0:1"));
    REQUIRE_THROWS(dur.parse("P0-1-0T0:0:0"));
    REQUIRE_THROWS(dur.parse(""));
    REQUIRE(dur.parse("P0-1").seconds() == 1 * days);
    REQUIRE(dur.parse("P0-0T0:0:1").seconds() == 1);
    REQUIRE(dur.parse("P-0-0T0:0:1").seconds() == -1);

    // arithmetic
    Duration yr("P1-0");
    Duration dy("P0-1");
    Duration s("P0-0T0:0:1");

    int daySeconds = 86400;
    int yearSeconds = 365 * daySeconds;
    REQUIRE((yr + dy).seconds() == yearSeconds + daySeconds);
    REQUIRE((dy - yr).seconds() == daySeconds - yearSeconds);

    yr += s;
    REQUIRE(yr.seconds() == yearSeconds + 1);
    yr -= s;
    REQUIRE(yr.seconds() == yearSeconds);
    yr *= 2;
    REQUIRE(yr.seconds() == 2 * yearSeconds);
    yr /= 2;
    REQUIRE(yr.seconds() == yearSeconds);

    // TimePoint manipulation
    TimePoint tt("2010-01-01T00:00:00Z");
    TimePoint tt_day("2010-01-02T00:00:00Z");

    REQUIRE(tt + dy == tt_day);
    REQUIRE(dy + tt == tt_day);

    TimePoint ttOther(tt);
    ttOther += dy;
    REQUIRE(ttOther == tt_day);

    Duration maybeADay = tt_day - tt;
    REQUIRE(maybeADay.seconds() == 1 * days);

    Duration tenSeconds(10.);
    REQUIRE(tenSeconds.seconds() == 10);

    Duration hundredSeconds("100.");
    REQUIRE(hundredSeconds.seconds() == 100);

    Duration minusThousandSeconds("-1000.");
    REQUIRE(minusThousandSeconds.seconds() == -1000);

    Duration myriadSeconds(1e4);
    REQUIRE(myriadSeconds.seconds() == 10000);

    Duration strMyriadSeconds("1e4");
    REQUIRE(strMyriadSeconds.seconds() == 10000);

    Duration complexStrSec("-1440.42e-2");
    REQUIRE(complexStrSec.seconds() == -14);

    Duration aDay(daySeconds);
    REQUIRE(aDay.seconds() == daySeconds);
}

TEST_CASE("gmtime and doy")
{
    TimePoint janfirst("2010-01-01T00:00:00Z");
    std::tm* timeStruct = janfirst.gmtime();
    REQUIRE(timeStruct->tm_yday == 0);

    // Test that leap years work
    TimePoint marchfirst("2010-03-01T00:00:00Z");
    REQUIRE(marchfirst.gmtime()->tm_yday == 59);
    TimePoint bissextile("2020-03-01T00:00:00Z");
    REQUIRE(bissextile.gmtime()->tm_yday == 60);
}

TEST_CASE("System time round trip")
{
    // Tests only system time library behaviour around leap days
    std::tm tm;
    int year = 1900;
    size_t month = 3;
    size_t day = 1;
    size_t hour = 0;
    size_t minute = 0;
    double second = 0.;

    std::stringstream(std::to_string(year) + "-" +
            std::to_string(month) + "-" +
            std::to_string(day) + "T" +
            std::to_string(hour) + ":" +
            std::to_string(minute) + ":" +
            std::to_string(second)) >> std::get_time(&tm, "%Y-%m-%dT%H:%M:%S");
    CHECK(tm.tm_year + 1900 == year);
    CHECK(tm.tm_mon == month-1);
    CHECK(tm.tm_mday == day);
    CHECK(tm.tm_hour == hour);
    CHECK(tm.tm_min == minute);
    CHECK(tm.tm_sec == second);

}

TEST_CASE("TimePoint date components")
{
    // A typical date
    testTimePointTime(2010, 1, 1, 0, 0, 0.);
    testTimePointTime(1899, 12, 31, 23, 59, 59);
    testTimePointTime(1900, 1, 1, 0, 0, 0);
    testTimePointTime(1900, 2, 28, 23, 59, 59);
    // Test a non-existent date outside the helper function
    TimePoint illegalDay("1900-02-29T00:00:00");
    REQUIRE(illegalDay.year() == 1900);
    REQUIRE(illegalDay.month() != 2);
    REQUIRE(illegalDay.day() != 29);
    REQUIRE(illegalDay.hour() == 0);
    REQUIRE(illegalDay.minute() == 0);
    REQUIRE(illegalDay.second() == 0);
    // We dare not venture further into 1900
    testTimePointTime(1901, 1, 1, 0, 0, 0);
    testTimePointTime(1901, 3, 1, 0, 0, 0);
    testTimePointTime(1904, 2, 28, 23, 59, 59);
    testTimePointTime(1904, 2, 29, 0, 0, 0);
    testTimePointTime(1904, 3, 1, 0, 0, 0);
    testTimePointTime(1969, 12, 31, 23, 59, 59);
    testTimePointTime(1970, 1, 1, 0, 0, 0);
    testTimePointTime(1999, 12, 31, 23, 59, 59);
    testTimePointTime(2000, 1, 1, 0, 0, 0);
    testTimePointTime(2000, 2, 28, 23, 59, 59);
    testTimePointTime(2000, 2, 29, 0, 0, 0);
    testTimePointTime(2000, 3, 1, 0, 0, 0);
}
TEST_SUITE_END();

}
