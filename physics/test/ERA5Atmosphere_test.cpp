/*!
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "include/ERA5Atmosphere.hpp"

#include "include/Time.hpp"

#include <string>

#include <ncDim.h>
#include <ncDouble.h>
#include <ncFile.h>
#include <ncVar.h>

namespace Nextsim {
// Signatures of the helper functions
using era5Buffer = Eigen::Array<double, Eigen::Dynamic, 1>;
std::string e5FilenameFromYear(const std::string& era5Name, size_t year);
era5Buffer getFileIndexData(const std::string& filename, size_t tIndex);
era5Buffer getVarIndexData(const std::string& era5Name, size_t year, size_t tIndex);
size_t timeIndexFromTM(const std::tm* tm);
era5Buffer getVarTimeData(const std::string& era5Name, const TimePoint& time);
ModelArray maFromERA5Buffer(const era5Buffer& buffer);
era5Buffer era5BufferHypot(const era5Buffer& x, const era5Buffer& y);

#ifndef TEST_FILES_DIR
#define TEST_FILES_DIR "."
#endif

const std::string testFilesDir = TEST_FILES_DIR;

int testFileYear = 2000;
const std::string era5NameLat = "lat";
const std::string latFileName = "ERA5_lat_y2000.nc";
const std::string era5NameLon = "lon";
const std::string lonFileName = "ERA5_lon_y2000.nc";
const std::string era5NameTimeTest = "timet";
const std::string timeTestFileName = "ERA5_timet_y2000.nc";

size_t nx = 1440;
size_t ny = 265;
double dlon = 0.25;
double lon0 = 0.0;
double dlat = -0.25;
double lat0 = 90.;
size_t nt = 8760;
double t0 = 0.;
double dt = 1.;

// Create longitude and latitude unprojected data files.
void createERA5TestFiles()
{
    // Spatial interpolation test files
    netCDF::NcFile latFile(testFilesDir + "/" + latFileName, netCDF::NcFile::replace, netCDF::NcFile::nc4);
    netCDF::NcFile lonFile(testFilesDir + "/" + lonFileName, netCDF::NcFile::replace, netCDF::NcFile::nc4);

    era5Buffer longitudeDim(nx, 1);
    era5Buffer latitudeDim(ny, 1);
    era5Buffer timeDim(nt, 1);
    era5Buffer lat2d(ny*nx, 1);
    era5Buffer lon2d(ny*nx, 1);

    for (size_t i = 0; i < nx; ++i) {
        longitudeDim(i, 0) = lon0 + dlon * i;
    }

    for (size_t j = 0; j < ny; ++j) {
        latitudeDim(j, 0) = lat0 + dlat * j;
        for (size_t i = 0; i < nx; ++i) {
            size_t idx = i + nx*(j);
            lon2d(idx, 0) = longitudeDim(i, 0);
            lat2d(idx, 0) = latitudeDim(j, 0);
        }
    }

    for (size_t t = 0; t < nt; ++t) {
        timeDim(t, 0) = t0 + t * dt;
    }

    std::vector<std::pair<netCDF::NcFile*, era5Buffer*>> fileDataPairs = {{&lonFile, &lon2d}, {&latFile, &lat2d}};
    for (std::pair<netCDF::NcFile*, era5Buffer*> fileDataPair : fileDataPairs) {
        netCDF::NcFile& file = *(fileDataPair.first);
        era5Buffer& data = *(fileDataPair.second);
        netCDF::NcDim lonDim = file.addDim("longitude", nx);
        netCDF::NcDim latDim = file.addDim("latitude", ny);
        netCDF::NcDim timDim = file.addDim("time", nt);
        file.addVar("longitude", netCDF::ncDouble, {lonDim,}).putVar({0,}, {nx,}, longitudeDim.data());
        file.addVar("latitude", netCDF::ncDouble, {latDim,}).putVar({0,}, {ny,}, latitudeDim.data());
        netCDF::NcVar timeVar(file.addVar("time", netCDF::ncDouble, {timDim,}));
        timeVar.putVar({0,}, {nt,}, timeDim.data());
        std::string timeUnits = "hours since 2010-01-01 00:00:00";
        timeVar.putAtt(std::string("units"), timeUnits);

        netCDF::NcVar dataVar(file.addVar("data", netCDF::ncDouble, {timDim, latDim, lonDim}));
        for (size_t t = 0; t < nt; ++t) {
            dataVar.putVar({t, 0, 0}, {1, ny, nx}, data.data());
        }
        file.close();
    }
}

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
TEST_CASE("hypot")
{
    era5Buffer a;
    era5Buffer b;
    era5Buffer c;
    a.resize(3*2, 1);
    b.resizeLike(a);
    c.resizeLike(a);
    a << 3, 5, 2.1, 1.5, 7, 5.5;
    b << 4, 12, 2, 0.8, 24, 4.8;
    c << 5, 13, 2.9, 1.7, 25, 7.3;

    era5Buffer d = era5BufferHypot(a, b);
    for (size_t i = 0; i < c.cols(); ++i) {
        REQUIRE(c(i,0) == doctest::Approx(d(i,0)).epsilon(1e-12));
    }
}

TEST_CASE("Interpolation tests")
{
    createERA5TestFiles();
}
}TEST_SUITE_END();
