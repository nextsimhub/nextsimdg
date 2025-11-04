/*!
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "include/TOPAZOcean.hpp"

#include "include/NetCDFForcings.hpp"
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
using Buffer = Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;

std::string topazFilenameFromYearMonth(const std::string& topazName, size_t year, size_t month);
ModelArray getTOPAZVarIndexData(const std::string& topazName, size_t year, size_t month, size_t day, const ModelArray& modelLon, const ModelArray& modelLat);

#ifndef TEST_FILES_DIR
#define TEST_FILES_DIR "."
#endif

const std::string testFilesDir = TEST_FILES_DIR;

size_t nx = 761;
size_t ny = 1101;
double dRad = radians(0.09);
double lon0 = 0.0;
double lat0 = radians(90.);
size_t nt = 1.;
double t0 = 438288; // hours between 1950-01-01T00:00:00 and 2000-01-01T00:00:00
double dt = 24.;

TEST_SUITE_BEGIN("TOPAZOcean");
TEST_CASE("topazFilename")
{
    std::string sst = "sst";
    std::string u = "u";

    size_t early = 1980;
    size_t recent = 2024;

    std::string testPath = "/test/path";
    TOPAZOcean::setDirectory(testPath);
    REQUIRE(topazFilenameFromYearMonth(sst, early, 1) == testPath + "/" + "TP4DAILY_198001_3m.nc");
    REQUIRE(topazFilenameFromYearMonth(sst, early, 12) == testPath + "/" + "TP4DAILY_198012_3m.nc");
    REQUIRE(topazFilenameFromYearMonth(u, recent, 1) == testPath + "/" + "TP4DAILY_202401_30m.nc");
}

TEST_CASE("Spatial interpolation from files")
{
    TOPAZOcean::setDirectory(testFilesDir);
    size_t fileYear = 2000;
    size_t fileMonth = 1;
    std::string fileDepth = "3m";
    std::string coordFileName = "TP4DAILY_200001_3m.nc";
    netCDF::NcFile coordFile(TOPAZOcean::addDirectory(coordFileName), netCDF::NcFile::replace, netCDF::NcFile::nc4);

    size_t nt = 1;
    Buffer timeDim(nt, 1);
    Buffer lon2d(nx, ny);
    Buffer lat2d(nx, ny);

    int ic = nx/2;
    int jc = ny/2;

    for (int j = 0; j < ny; ++j) {
        double y = (j - jc) * dRad;
        for (int i = 0; i < nx; ++i) {
            double x = (i - ic) * dRad;
            double rho = sqrt(x*x + y*y);
            double c = 2 * atan(rho / 2);
            lat2d(i, j) = degrees(asin(cos(c)));
            lon2d(i, j) = degrees(atan2(x*sin(c), rho*cos(lat0)*cos(c) - y*sin(lat0)*sin(c)));
        }
    }

    for (size_t t = 0; t < nt; ++t) {
        timeDim(t, 0) = t0 + t * dt;
    }

    netCDF::NcFile& file = coordFile;
    netCDF::NcDim xDim = file.addDim("x", nx);
    netCDF::NcDim yDim = file.addDim("y", ny);
    netCDF::NcDim timDim = file.addDim("time", nt);
    netCDF::NcVar lonCoordVar(file.addVar("longitude", netCDF::ncDouble, {yDim, xDim}));
    lonCoordVar.putVar({0, 0}, {ny, nx}, lon2d.data());
    netCDF::NcVar latCoordVar(file.addVar("latitude", netCDF::ncDouble, {yDim, xDim}));
    latCoordVar.putVar({0, 0}, {ny, nx}, lat2d.data());
    netCDF::NcVar timeVar(file.addVar("time", netCDF::ncDouble, {timDim,}));
    timeVar.putVar({0,}, {nt,}, timeDim.data());
    // Longitude and latitude as data, rather than coordinates
    netCDF::NcVar lonVar(file.addVar("lon", netCDF::ncDouble, {timDim, yDim, xDim}));
    lonVar.putVar({0, 0, 0}, {1, ny, nx}, lon2d.data());
    netCDF::NcVar latVar(file.addVar("lat", netCDF::ncDouble, {timDim, yDim, xDim}));
    latVar.putVar({0, 0, 0}, {1, ny, nx}, lat2d.data());
    std::string timeUnits = "hours since 1950-01-01 00:00:00";
    timeVar.putAtt(std::string("units"), timeUnits);
    timeVar.putAtt(std::string("calendar"), std::string("standard"));

    file.close();

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
        ModelArray::setDimensions(ModelArray::Type::H, {nxt, nyt});
        ModelArray testLon = getTOPAZVarIndexData("lon", fileYear, fileMonth, 1, lonTarg, latTarg);
        ModelArray testLat = getTOPAZVarIndexData("lat", fileYear, fileMonth, 1, lonTarg, latTarg);
        std::vector<size_t> test_i = {0, 1, 95, 152, 153};
        std::vector<size_t> test_j = {0, 1, 36, 119, 120};

        double prec = 1;
        for (size_t j : test_j) {
            for (size_t i : test_i) {
                REQUIRE(std::fmod(testLon(i, j) - lonTarg(i, j), 360.) == doctest::Approx(0.).epsilon(prec));
                REQUIRE(testLat(i, j) == doctest::Approx(latTarg(i, j)).epsilon(prec));
            }
        }
    }
    if (std::filesystem::exists(TOPAZOcean::addDirectory(coordFileName))) {
        std::filesystem::remove(TOPAZOcean::addDirectory(coordFileName));
    }
}
TEST_SUITE_END();
}
