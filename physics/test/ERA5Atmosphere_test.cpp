/*!
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "include/ERA5Atmosphere.hpp"

#include "include/Time.hpp"

#include <string>

namespace Nextsim {
// Signatures of the helper functions
using era5Buffer = Eigen::Array<double, Eigen::Dynamic, 1>;
std::string e5FilenameFromYear(const std::string& era5Name, size_t year);
era5Buffer getFileIndexData(const std::string& filename, size_t tIndex);
era5Buffer getVarIndexData(const std::string& era5Name, size_t year, size_t tIndex);
size_t timeIndexFromTM(const std::tm* tm);
era5Buffer getVarTimeData(const std::string& era5Name, const TimePoint& time);
ModelArray maFromERA5Buffer(const era5Buffer& buffer);

TEST_SUITE_BEGIN("ERA5Atmosphere");
TEST_CASE("e5Filename")
{
    std::string t2m = "t2m";
    std::string downwardLongwaveFlux = "msdwlwrf";

    size_t early = 1980;
    size_t recent = 2024;

    REQUIRE(e5FilenameFromYear(t2m, early) == "ERA5_t2m_y1980.nc");
    REQUIRE(e5FilenameFromYear(downwardLongwaveFlux, recent) == "ERA5_msdwlwrf_y2024.nc");
}
TEST_CASE("Time index")
{
    std::tm tm1{};
    // 00:00 1st January 2010
    tm1.tm_year = 2010 - 1900;
    tm1.tm_mon = 0; // Month is 0 based
    tm1.tm_yday = 0; // Day of the month is 1-based
    tm1.tm_hour = 0;
    tm1.tm_min = 0;
    tm1.tm_sec = 0;

    REQUIRE(timeIndexFromTM(&tm1) == 0);
    ++tm1.tm_sec;
    REQUIRE(timeIndexFromTM(&tm1) == 0);
    ++tm1.tm_min;
    --tm1.tm_sec;
    REQUIRE(timeIndexFromTM(&tm1) == 0);
    ++tm1.tm_hour;
    --tm1.tm_min;
    REQUIRE(timeIndexFromTM(&tm1) == 1);
    ++tm1.tm_yday;
    --tm1.tm_hour;
    REQUIRE(timeIndexFromTM(&tm1) == 24);
    --tm1.tm_yday;
    ++tm1.tm_year;
    REQUIRE(timeIndexFromTM(&tm1) == 0);
}
}
#undef ERA5TESTING
