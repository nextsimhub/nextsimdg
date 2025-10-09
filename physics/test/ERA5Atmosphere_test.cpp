/*!
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "include/ERA5Atmosphere.hpp"

#include "include/ModelArray.hpp"
#include "include/Time.hpp"

#include <filesystem>
#include <cmath>
#include <string>

#include <ncDim.h>
#include <ncDouble.h>
#include <ncFile.h>
#include <ncVar.h>

namespace Nextsim {
// Signatures of the helper functions
using era5Buffer = Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
std::string e5FilenameFromYear(const std::string& era5Name, size_t year);
era5Buffer getFileIndexData(const std::string& filename, size_t tIndex);
era5Buffer getVarIndexData(const std::string& era5Name, size_t year, size_t tIndex);
size_t timeIndexFromTM(const std::tm* tm);
era5Buffer getVarTimeData(const std::string& era5Name, const TimePoint& time);
ModelArray maFromERA5Buffer(const era5Buffer& buffer, const ModelArray& destLon, const ModelArray& destLat);
era5Buffer era5BufferHypot(const era5Buffer& x, const era5Buffer& y);

#ifndef TEST_FILES_DIR
#define TEST_FILES_DIR "."
#endif

const std::string testFilesDir = TEST_FILES_DIR;

int testFileYear = 2000;
const std::string era5NameTime = "timet";
const std::string timeFileName = "ERA5_timet_y2000.nc";

size_t nx = 1440;
size_t ny = 265;
double dlon = 0.25;
double lon0 = 0.0;
double dlat = -0.25;
double lat0 = 90.;
size_t nt = 24;
double t0 = 0.;
double dt = 1.;

static const double pi = 0x3.243F68885Ap0;
double radians(double deg)
{
    return deg * pi / 180;
}
double degrees(double rad)
{
    return rad * 180. / pi;
}
// Create longitude and latitude unprojected data files.
void createERA5TestFiles()
{
    // Temporal interpolation test file
    netCDF::NcFile timeFile(testFilesDir + "/" + timeFileName, netCDF::NcFile::replace, netCDF::NcFile::nc4);

    era5Buffer longitudeDim(nx, 1);
    era5Buffer latitudeDim(ny, 1);
    era5Buffer timeDim(nt, 1);
    era5Buffer time2d(nx, ny);

    for (size_t i = 0; i < nx; ++i) {
        longitudeDim(i, 0) = lon0 + dlon * i;
    }

    for (size_t j = 0; j < ny; ++j) {
        latitudeDim(j, 0) = lat0 + dlat * j;
//        }
    }

    for (size_t t = 0; t < nt; ++t) {
        timeDim(t, 0) = t0 + t * dt;
    }

    time2d.setZero();
    netCDF::NcDim lonDim = timeFile.addDim("longitude", nx);
    netCDF::NcDim latDim = timeFile.addDim("latitude", ny);
    netCDF::NcDim timDim = timeFile.addDim("time", nt);
    timeFile.addVar("longitude", netCDF::ncDouble, {lonDim,}).putVar({0,}, {nx,}, longitudeDim.data());
    timeFile.addVar("latitude", netCDF::ncDouble, {latDim,}).putVar({0,}, {ny,}, latitudeDim.data());
    netCDF::NcVar timeVar(timeFile.addVar("time", netCDF::ncDouble, {timDim,}));
    timeVar.putVar({0,}, {nt,}, timeDim.data());
    std::string timeUnits = "hours since 1900-01-01 00:00:00";
    timeVar.putAtt(std::string("units"), timeUnits);
    timeVar.putAtt(std::string("calendar"), std::string("standard"));

    netCDF::NcVar dataVar(timeFile.addVar("data", netCDF::ncDouble, {timDim, latDim, lonDim}));
    for (size_t t = 0; t < nt; ++t) {
        dataVar.putVar({t, 0, 0}, {1, ny, nx}, time2d.data());
        time2d += 1;
    }
    timeFile.close();
}

TEST_SUITE_BEGIN("ERA5Atmosphere");
TEST_CASE("e5Filename")
{
    std::string t2m = "t2m";
    std::string downwardLongwaveFlux = "msdwlwrf";

    size_t early = 1980;
    size_t recent = 2024;

    std::string testPath = "/test/path";
    ERA5Atmosphere::setDirectory(testPath);
    REQUIRE(e5FilenameFromYear(t2m, early) == testPath + "/" + "ERA5_t2m_y1980.nc");
    REQUIRE(e5FilenameFromYear(downwardLongwaveFlux, recent) == testPath + "/" + "ERA5_msdwlwrf_y2024.nc");
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

TEST_CASE("Interpolation tests I: time index from file")
{
    ERA5Atmosphere::setDirectory(testFilesDir);
    if (!std::filesystem::exists(ERA5Atmosphere::addDirectory(timeFileName))) {
        createERA5TestFiles();
    }

    const size_t tIdx = 4;
    era5Buffer timeData = getFileIndexData(ERA5Atmosphere::addDirectory(timeFileName), tIdx);
    REQUIRE(timeData(0, 0) == double(tIdx));
}

TEST_CASE("Interpolation tests II: variable data from year and index")
{
    if (!std::filesystem::exists(testFilesDir + "/" + timeFileName)) {
        createERA5TestFiles();
    }
    ERA5Atmosphere::setDirectory(testFilesDir);
    REQUIRE(ERA5Atmosphere::getDirectory() == testFilesDir);

    size_t timeIdx = 5;
    era5Buffer timeData = getVarIndexData(era5NameTime, 2000, timeIdx);
    REQUIRE(timeData(0, 0) == double(timeIdx));

    // No timet file should exist for 1999
    timeIdx = 7;
    REQUIRE_THROWS(timeData = getVarIndexData(era5NameTime, 1999, timeIdx));

}

TEST_CASE("Interpolation tests III: variable data from TimePoint")
{
    if (!std::filesystem::exists(testFilesDir + "/" + timeFileName)) {
        createERA5TestFiles();
    }
    ERA5Atmosphere::setDirectory(testFilesDir);
    // A quarter of the way between index 6 and index 7
    TimePoint qpt("2000-01-01T06:15:00Z");
    double targetValue = 6.25;
    era5Buffer timeData = getVarTimeData(era5NameTime, qpt);
    REQUIRE(timeData(0, 0) == targetValue);
}

TEST_CASE("Interpolation tests IV: spatial interpolation")
{
    using std::sqrt;
    using std::cos;
    using std::sin;
    using std::atan;
    using std::atan2;
    using std::asin;

    era5Buffer lon(nx, 1);
    era5Buffer lat(ny, 1);
    era5Buffer lat2d(nx, ny);
    era5Buffer lon2d(nx, ny);

    for (size_t i = 0; i < nx; ++i) {
        lon(i, 0) = lon0 + dlon * i;
    }

    for (size_t j = 0; j < ny; ++j) {
        lat(j, 0) = lat0 + dlat * j;
        for (size_t i = 0; i < nx; ++i) {
            lon2d(i, j) = lon(i, 0);
            lat2d(i, j) = lat(j, 0);
        }
    }

    // Longitude and latitude of the target grid
    size_t nxt = 154;
    size_t nyt = 121;
    double lonCentreDeg = 180.;
    double lonC = radians(lonCentreDeg);
    double latCentreDeg = 82.;
    double latC = radians(latCentreDeg);
    double dDeg = 0.25; // Resolution in degrees

    ModelArray::setDimension(ModelArray::Dimension::X, nxt);
    ModelArray::setDimension(ModelArray::Dimension::Y, nyt);
    ModelArray lonTarg(ModelArray::Type::H);
    ModelArray latTarg(ModelArray::Type::H);
    lonTarg.resize();
    latTarg.resize();
    // Create the target longitude & latitude arrays: polar stereographic
    int ic = nxt/2;
    int jc = nyt/2;

    for (int j = 0; j < nyt; ++j) {
        double y = (j - jc) * radians(dDeg);
        for (int i = 0; i < nxt; ++i) {
            double x = (i - ic) * radians(dDeg);
            double rho = sqrt(x*x + y*y);
            double c = 2 * atan(rho / 2);
            latTarg(i, j) = degrees(asin(cos(c)*sin(latC) + y*sin(c)*cos(latC)/rho));
            lonTarg(i, j) = degrees(lonC + atan2(x*sin(c), rho*cos(latC)*cos(c) - y*sin(latC)*sin(c)));
        }
    }
    ModelArray testLat(maFromERA5Buffer(lat2d, lonTarg, latTarg));
    size_t testi = 20;
    size_t testj = 45;
    REQUIRE(testLat(testi, testj) == latTarg(testi, testj));

    ModelArray testLon(maFromERA5Buffer(lon2d, lonTarg, latTarg));
    testi = 45;
    testj = 20;
    REQUIRE(testLon(testi, testj) == lonTarg(testi, testj));
}
TEST_SUITE_END();
}
