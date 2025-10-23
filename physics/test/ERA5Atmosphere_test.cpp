/*!
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "include/ERA5Atmosphere.hpp"

#include "include/ModelArray.hpp"
#include "include/Time.hpp"

#include <cmath>
#include <cstdint>
#include <filesystem>
#include <string>
#include <tuple>
#include <vector>

#include <ncDim.h>
#include <ncDouble.h>
#include <ncFile.h>
#include <ncShort.h>
#include <ncVar.h>

namespace Nextsim {
// Signatures of the helper functions
using era5Buffer = Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;
using era5ShortBuffer = Eigen::Array<std::int16_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;

std::string e5FilenameFromYear(const std::string& era5Name, size_t year);
era5Buffer getFileIndexData(const std::string& filename, size_t tIndex);
era5Buffer getVarIndexData(const std::string& era5Name, size_t year, size_t tIndex);
size_t timeIndex(const TimePoint& tm);
era5Buffer getVarTimeData(const std::string& era5Name, const TimePoint& time);
ModelArray maFromERA5Buffer(const era5Buffer& buffer, const ModelArray& destLon, const ModelArray& destLat);
era5Buffer era5BufferHypot(const era5Buffer& x, const era5Buffer& y);
std::pair<double, double> getLatitudeCoeffs(const std::string& filename);

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
void createERA5TimeTestFiles()
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
    TimePoint time("2010-01-01T00:00:00");

    REQUIRE(timeIndex(time) == 0);
    REQUIRE(timeIndex(time + Duration(1)) == 0);
    REQUIRE(timeIndex(time + Duration(60)) == 0);
    REQUIRE(timeIndex(time + Duration(3600)) == 1);
    REQUIRE(timeIndex(time + Duration(86400)) == 24);
    REQUIRE(timeIndex(time + Duration("P0-1")) == 24);
    REQUIRE(timeIndex(time + Duration("P1-0")) == 0);
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

TEST_CASE("Time interpolation tests")
{
    ERA5Atmosphere::setDirectory(testFilesDir);
    if (!std::filesystem::exists(ERA5Atmosphere::addDirectory(timeFileName))) {
        createERA5TimeTestFiles();
    }

    SUBCASE("I: time index from file") {
        const size_t tIdx = 4;
        era5Buffer timeData = getFileIndexData(ERA5Atmosphere::addDirectory(timeFileName), tIdx);
        REQUIRE(timeData(0, 0) == double(tIdx));
    }

    SUBCASE("II: variable data from year and index")
    {
        if (!std::filesystem::exists(testFilesDir + "/" + timeFileName)) {
        createERA5TimeTestFiles();
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

    SUBCASE("III: variable data from TimePoint")
    {
        if (!std::filesystem::exists(testFilesDir + "/" + timeFileName)) {
            createERA5TimeTestFiles();
        }
        ERA5Atmosphere::setDirectory(testFilesDir);
        // A quarter of the way between index 6 and index 7
        TimePoint qpt("2000-01-01T06:15:00Z");
        double targetValue = 6.25;
        era5Buffer timeData = getVarTimeData(era5NameTime, qpt);
        REQUIRE(timeData(0, 0) == targetValue);
    }

    if (std::filesystem::exists(testFilesDir + "/" + timeFileName)) {
        std::filesystem::remove(testFilesDir + "/" + timeFileName);
    }
}
TEST_CASE("Spatial interpolation from field")
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

    double prec = 1e-14;
    ModelArray testLat(maFromERA5Buffer(lat2d, lonTarg, latTarg));
    ModelArray testLon(maFromERA5Buffer(lon2d, lonTarg, latTarg));
    size_t testi = 20;
    size_t testj = 45;
    REQUIRE(testLat(testi, testj) == doctest::Approx(latTarg(testi, testj)).epsilon(prec));
    REQUIRE(testLon(testi, testj) == doctest::Approx(lonTarg(testi, testj)).epsilon(prec));

    testi = 45;
    testj = 20;
    REQUIRE(testLat(testi, testj) == doctest::Approx(latTarg(testi, testj)).epsilon(prec));
    REQUIRE(testLon(testi, testj) == doctest::Approx(lonTarg(testi, testj)).epsilon(prec));
}

TEST_CASE("Spatial interpolation from files")
{
    ERA5Atmosphere::setDirectory(testFilesDir);
    std::string lonFileName = "ERA5_lonx_y2000.nc";
    netCDF::NcFile lonFile(ERA5Atmosphere::addDirectory(lonFileName), netCDF::NcFile::replace, netCDF::NcFile::nc4);
    std::string latFileName = "ERA5_laty_y2000.nc";
    netCDF::NcFile latFile(ERA5Atmosphere::addDirectory(latFileName), netCDF::NcFile::replace, netCDF::NcFile::nc4);

    size_t nt = 1;
    era5Buffer longitudeDim(nx, 1);
    era5Buffer latitudeDim(ny, 1);
    era5Buffer timeDim(nt, 1);
    era5ShortBuffer lon2d(nx, ny);
    era5ShortBuffer lat2d(nx, ny);

    double aLon = 0.25;
    double bLon = 0.;
    double aLat = -0.01;
    double bLat = 90.;

    for (size_t i = 0; i < nx; ++i) {
        double lon = lon0 + dlon * i;
        longitudeDim(i, 0) = lon;
        lon2d(i, Eigen::all) = static_cast<std::int16_t>((lon - bLon) / aLon);
    }

    for (size_t j = 0; j < ny; ++j) {
        double lat = lat0 + dlat * j;
        latitudeDim(j, 0) = lat;
        lat2d(Eigen::all, j) = static_cast<std::int16_t>((lat - bLat) / aLat);
    }

    for (size_t t = 0; t < nt; ++t) {
        timeDim(t, 0) = t0 + t * dt;
    }

    for (std::tuple<netCDF::NcFile*, era5ShortBuffer*, double, double> fileData : std::vector<std::tuple<netCDF::NcFile*, era5ShortBuffer*, double, double>> {
            {&lonFile, &lon2d, aLon, bLon},
            {&latFile, &lat2d, aLat, bLat},
        }) {
        netCDF::NcFile& file = *(std::get<0>(fileData));
        netCDF::NcDim lonDim = file.addDim("longitude", nx);
        netCDF::NcDim latDim = file.addDim("latitude", ny);
        netCDF::NcDim timDim = file.addDim("time", nt);
        file.addVar("longitude", netCDF::ncDouble, {lonDim,}).putVar({0,}, {nx,}, longitudeDim.data());
        file.addVar("latitude", netCDF::ncDouble, {latDim,}).putVar({0,}, {ny,}, latitudeDim.data());
        netCDF::NcVar timeVar(file.addVar("time", netCDF::ncDouble, {timDim,}));
        timeVar.putVar({0,}, {nt,}, timeDim.data());
        std::string timeUnits = "hours since 1900-01-01 00:00:00";
        timeVar.putAtt(std::string("units"), timeUnits);
        timeVar.putAtt(std::string("calendar"), std::string("standard"));

        netCDF::NcVar dataVar(file.addVar("data", netCDF::ncShort, {timDim, latDim, lonDim}));
        dataVar.putAtt(std::string("scale_factor"), netCDF::ncDouble, std::get<2>(fileData));
        dataVar.putAtt(std::string("add_offset"), netCDF::ncDouble, std::get<3>(fileData));
        for (size_t t = 0; t < nt; ++t) {
            dataVar.putVar({t, 0, 0}, {1, ny, nx}, std::get<1>(fileData)->data());
        }
        file.close();
    }

    SUBCASE("Lon and lat dimensions") {
        auto [dlat, lat0] = getLatitudeCoeffs(ERA5Atmosphere::addDirectory(lonFileName));
        REQUIRE(lat0 == 90.);
        REQUIRE(dlat == -0.25);
    }
    SUBCASE("Interpolation") {
        // Create the coordinates of the target grid
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
        era5Buffer readLon(getFileIndexData(ERA5Atmosphere::addDirectory("ERA5_lonx_y2000.nc"), 0));
        era5Buffer readLat(getFileIndexData(ERA5Atmosphere::addDirectory("ERA5_laty_y2000.nc"), 0));
        ModelArray::setDimensions(ModelArray::Type::H, {nxt, nyt});
        ModelArray testLon;
        testLon.resize();
        testLon = maFromERA5Buffer(readLon, lonTarg, latTarg);
        ModelArray testLat;
        testLat.resize();
        testLat = maFromERA5Buffer(readLat, lonTarg, latTarg);
        std::vector<size_t> test_i = {0, 1, 95, 152, 153};
        std::vector<size_t> test_j = {0, 1, 36, 119, 120};

        double prec = 1e-14;
        for (size_t j : test_j) {
            for (size_t i : test_i) {
                REQUIRE(testLon(i, j) == doctest::Approx(lonTarg(i, j)).epsilon(prec));
                REQUIRE(testLat(i, j) == doctest::Approx(latTarg(i, j)).epsilon(prec));
            }
        }
    }
    if (std::filesystem::exists(ERA5Atmosphere::addDirectory(lonFileName))) {
        std::filesystem::remove(ERA5Atmosphere::addDirectory(lonFileName));
    }
    if (std::filesystem::exists(ERA5Atmosphere::addDirectory(latFileName))) {
        std::filesystem::remove(ERA5Atmosphere::addDirectory(latFileName));
    }

}
TEST_SUITE_END();
}
