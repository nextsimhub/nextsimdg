/*!
 * @file Time_test.cpp
 *
 * @date Jun 7, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/Time.hpp"

#include <sstream>

#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

namespace Nextsim {

TEST_CASE("TimePoint parsing and formating", "[TimePoint]")
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

TEST_CASE("Julian-Gregorian shifts", "[Time]")
{
    // In the 1500s, the inital 10 day shift.
    REQUIRE(julianGregorianShiftDays(1582) == -10);
    // In the 1600s, 1600 had been a leap year, as expected, so still 10 days.
    REQUIRE(julianGregorianShiftDays(1688) == -10);
    // In 1700, the different leap day had not yet been reached.
    REQUIRE(julianGregorianShiftDays(1700) == -10);
    // By 1701, there had been one extra skipped Gregorian leap day.
    REQUIRE(julianGregorianShiftDays(1701) == -11);
    // By 1970, there had been 3 further fewer Gregorian leap days (1700, 1800, 1900).
    REQUIRE(julianGregorianShiftDays(1970) == -13);
    // In 2525, if man is still alive, there will have been 7 further fewer
    // Gregorian leap days (1700, 1800, 1900, 2100, 2200, 2300, 2500).
    REQUIRE(julianGregorianShiftDays(2525) == -17);

}

TEST_CASE("mkgmtime", "[Time]")
{
    const int days = 24 * 60 * 60;

    const char* iso = TimePoint::ymdhmsFormat.c_str();

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

TEST_CASE("timeFromISO", "[Time]")
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

TEST_CASE("TimePoints", "[TimePoint]")
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

TEST_CASE("Durations", "[Duration]")
{
    const int days = 24 * 60 * 60;

    Duration dur;
    // Basic values
    REQUIRE_THROWS(dur.parse("0-0-0T0:0:1"));
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
}
}
