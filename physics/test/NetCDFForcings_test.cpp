/*!
 * @file NetCDFForcings_test.cpp
 *
 * @date Oct 21, 2025
 * @author Tim Spain <timothy.spain@nersc.no>
 */
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "include/NetCDFForcings.hpp"

#include "include/ModelArray.hpp"

namespace Nextsim {
static const double pi = 0x3.243F68885Ap0;
double radians(double deg)
{
    return deg * pi / 180;
}
double degrees(double rad)
{
    return rad * 180. / pi;
}

TEST_CASE("Spatial interpolation from field")
{
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

    SUBCASE("ERA5-like test")
    {
        using std::sqrt;
        using std::cos;
        using std::sin;
        using std::atan;
        using std::atan2;
        using std::asin;

        size_t nx = 1440;
        size_t ny = 265;
        double dlon = 0.25;
        double lon0 = 0.0;
        double dlat = -0.25;
        double lat0 = 90.;

        NetCDFForcings::Buffer lon(nx, 1);
        NetCDFForcings::Buffer lat(ny, 1);
        NetCDFForcings::Buffer lat2d(nx, ny);
        NetCDFForcings::Buffer lon2d(nx, ny);

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

        ModelArray iFrac = lonTarg / dlon;
        ModelArray jFrac = (latTarg - lat0) / dlat;

        double prec = 1e-14;
        ModelArray testLat(NetCDFForcings::maFromBuffer(lat2d, iFrac, jFrac));
        ModelArray testLon(NetCDFForcings::maFromBuffer(lon2d, iFrac, jFrac));
        size_t testi = 20;
        size_t testj = 45;
        REQUIRE(testLat(testi, testj) == doctest::Approx(latTarg(testi, testj)).epsilon(prec));
        REQUIRE(testLon(testi, testj) == doctest::Approx(lonTarg(testi, testj)).epsilon(prec));

        testi = 45;
        testj = 20;
        REQUIRE(testLat(testi, testj) == doctest::Approx(latTarg(testi, testj)).epsilon(prec));
        REQUIRE(testLon(testi, testj) == doctest::Approx(lonTarg(testi, testj)).epsilon(prec));
    }

    SUBCASE("TOPAZ4-like test")
    {
        using std::sqrt;
        using std::cos;
        using std::sin;
        using std::atan;
        using std::atan2;
        using std::asin;

        // Centred polar stereographic project, 0.09˚ grid spacing (= 10km)
        size_t nx = 761;
        size_t ny = 1101;

        double dDeg = 0.09;
        double dRad = radians(dDeg);
        double lonCentreDeg = 0.;
        double lonC = radians(lonCentreDeg);
        double latCentreDeg = 90.;
        double latC = radians(latCentreDeg);

        NetCDFForcings::Buffer plon2d(nx, ny);
        NetCDFForcings::Buffer plat2d(nx, ny);
        // Create the target longitude & latitude arrays: polar stereographic
        int ic = nx/2;
        int jc = ny/2;

        for (int j = 0; j < ny; ++j) {
            double y = (j - jc) * dRad;
            for (int i = 0; i < nx; ++i) {
                double x = (i - ic) * dRad;
                double rho = sqrt(x*x + y*y);
                double c = 2 * atan(rho / 2);
                plat2d(i, j) = degrees(asin(cos(c)*sin(latC) + y*sin(c)*cos(latC)/rho));
                plon2d(i, j) = degrees(lonC + atan2(x*sin(c), rho*cos(latC)*cos(c) - y*sin(latC)*sin(c)));
            }
        }

        ModelArray iFrac;
        ModelArray jFrac;
        iFrac.resize();
        jFrac.resize();

        // Create the arrays of indices into these source arrays for the positions in the model arrays
        for (int j = 0; j < ModelArray::size(ModelArray::Dimension::Y); ++j) {
            for (int i = 0; i < ModelArray::size(ModelArray::Dimension::X); ++i) {
                double radLat = radians(latTarg(i, j));
                double rho = 2 * cos(radLat) / (1 + sin(radLat));
                double radLon = radians(lonTarg(i, j));
                iFrac(i, j) = rho * sin(radLon) / dRad + ic;
                jFrac(i, j) = -rho * cos(radLon) / dRad + jc;
            }
        }

        double prec = 1e-6;
        ModelArray testLat(NetCDFForcings::maFromBuffer(plat2d, iFrac, jFrac));
        ModelArray testLon(NetCDFForcings::maFromBuffer(plon2d, iFrac, jFrac));
        size_t testi = 20;
        size_t testj = 45;
        REQUIRE(testLat(testi, testj) == doctest::Approx(latTarg(testi, testj)).epsilon(prec));
        REQUIRE(testLon(testi, testj) == doctest::Approx(lonTarg(testi, testj)).epsilon(prec));

        testi = 45;
        testj = 20;
        REQUIRE(testLat(testi, testj) == doctest::Approx(latTarg(testi, testj)).epsilon(prec));
        REQUIRE(testLon(testi, testj) == doctest::Approx(lonTarg(testi, testj)).epsilon(prec));
    }

}

} /* namespace Nextsim */
