/*!
 * @file FieldAdvection_test.cpp
 *
 * @date Aug 1, 2024
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "include/AdvectionDynamicsKernel.hpp"
#include "include/ModelArray.hpp"
#include "include/Time.hpp"

#include <fstream>

namespace Nextsim {

void writeDataPGM(const ModelArray& data, const std::string& fileName, size_t k)
{
    std::ofstream pgm(fileName);

    // Magic number
    pgm << "P2" << std::endl;
    pgm << "# Ice temperature in negative centi-Celsius" << std::endl;
    size_t nx = ModelArray::size(ModelArray::Dimension::X);
    size_t ny = ModelArray::size(ModelArray::Dimension::Y);
    pgm << nx << " " << ny << std::endl;
    long maxVal = 255L;
    pgm << maxVal << std::endl;
    for (size_t j = 0; j < ny; ++j) {
        for (size_t i = 0; i < nx; ++i) {
            pgm << std::max(0L, std::min(maxVal, std::lround(data(i, j, k)))) << " ";
        }
        pgm << std::endl;
    }
    pgm.close();

}

std::pair<double, double> pseudoFT(const ModelArray& data, size_t x0, size_t y0, size_t r)
{
    // Sample 8 points around the circle to do a basic Fourier transform
    const size_t nSamples = 8;
    std::array<size_t, nSamples> xTest = { x0 + r, x0 + r, x0, x0 - r, x0 - r, x0 - r, x0, x0 + r };
    std::array<size_t, nSamples> yTest = { y0, y0 + r, y0 + r, y0 + r, y0, y0 - r, y0 - r, y0 - r };
    std::array<double, nSamples> cosFactor = { 1, 0, -1, 0, 1, 0, -1, 0 };
    std::array<double, nSamples> sinFactor = { 0, 1, 0, -1, 0, 1, 0, -1 };

    double cosTerm = 0;
    double sinTerm = 0;
    for (size_t theta = 0; theta < nSamples; ++theta) {
        cosTerm += cosFactor[theta] * data(xTest[theta], yTest[theta]);
        sinTerm += sinFactor[theta] * data(xTest[theta], yTest[theta]);
    }
    return {cosTerm, sinTerm};
}

void writeDataPGM(const ModelArray& data, const std::string& fileName)
{
    std::ofstream pgm(fileName);

    // Magic number
    pgm << "P2" << std::endl;
    pgm << "# Depth of snow in cm" << std::endl;
    size_t nx = ModelArray::size(ModelArray::Dimension::X);
    size_t ny = ModelArray::size(ModelArray::Dimension::Y);
    pgm << nx << " " << ny << std::endl;
    long maxVal = 255L;
    pgm << maxVal << std::endl;
    for (size_t j = 0; j < ny; ++j) {
        for (size_t i = 0; i < nx; ++i) {
            pgm << std::max(0L, std::min(maxVal, std::lround(data(i, j)))) << " ";
        }
        pgm << std::endl;
    }
    pgm.close();

}


TEST_SUITE_BEGIN("Field advection");
TEST_CASE("Advect a field")
{
    // Parameters of the mesh
    size_t nx = 101;
    size_t ny = 101;
    size_t nz = 3;

    double dx = 1000; //m
    double dy = 1000; //m
    double r0 = 40000; //m
    double deltaT = 25; //s

    // Create the grid coordinates
    ModelArray::setDimension(ModelArray::Dimension::X, nx);
    ModelArray::setDimension(ModelArray::Dimension::Y, ny);
    ModelArray::setDimension(ModelArray::Dimension::Z, nz);
    ModelArray::setDimension(ModelArray::Dimension::XVERTEX, nx + 1);
    ModelArray::setDimension(ModelArray::Dimension::YVERTEX, ny + 1);
    VertexField coords(ModelArray::Type::VERTEX);
    coords.resize();
    for (size_t j = 0; j <= ny; ++j) {
        for (size_t i = 0; i <= nx; ++i) {
            coords.components({i, j})[0] = dx * i;
            coords.components({i, j})[1] = dy * j;
        }
    }
    HField mask(ModelArray::Type::H);
    mask.resize();
    mask = 1.; // water/ice everywhere

    AdvectionDynamicsKernel<6> advection;
    advection.initialise(coords, false, mask);

    UField u(ModelArray::Type::U);
    VField v(ModelArray::Type::V);
    u.resize();
    v.resize();
    HField hice(ModelArray::Type::H);
    HField cice(ModelArray::Type::H);
    hice.resize();
    cice.resize();
    HField hsnow(ModelArray::Type::H);
    hsnow.resize();
    ZField tice(ModelArray::Type::Z);
    tice.resize();
    
    HField x(ModelArray::Type::H);
    HField y(ModelArray::Type::H);
    x.resize();
    y.resize();
    size_t i0 = nx / 2;
    size_t j0 = ny / 2;
    double x0 = (nx / 2.) * dx;
    double y0 = (ny / 2.) * dy;
    double omega = 6.28e-5;
    double hice0 = 1.;
    double cice0 = 1.;
    double hsnowA = 1.;
//    double tice0 = -1.;
    double ticeA = 1.;
    // Slope characteristics
    double rs = 1 * dx;
    for (size_t j = 0; j < ny; ++j) {
        for (size_t i = 0; i < nx; ++i) {
            // Linear coordinates
            double xij = i * dx - x0;
            double yij = j * dy - y0;
            x(i, j) = xij;
            y(i, j) = yij;
            // Square coordinates
            double xx = xij * xij;
            double yy = yij * yij;
            // Radial coordinates
            double rr = xx + yy;
            double r = std::sqrt(rr);
            // Hsnow function: cos^2 + 1
            double hsnow0 = hsnowA * xx / rr + 1;
            double tice0 = ticeA * xx / rr - 1;
            if (r < r0) {
                hice(i, j) = hice0;
                cice(i, j) = cice0;
                u(i, j) = -omega * yij;
                v(i, j) = omega * xij;
                hsnow(i, j) = hsnow0;
                for (size_t k = 0; k < nz; ++k) {
                    tice(i, j, k) = tice0;
                }
            } else {
                double slope = std::max((r0 + rs - r) / rs, 0.);
                hice(i, j) = hice0 * slope;
                cice(i, j) = cice0 * slope;
                // u, v are defined if there is any ice
                if (slope > 0.) {
                    u(i, j) = -omega * yij;
                    v(i, j) = omega * xij;
                } else {
                    u(i, j) = 0.;
                    v(i, j) = 0.;
                }
                hsnow(i, j) = hsnow0 * slope;
                for (size_t k = 0; k < nz; ++k) {
                    if (cice(i, j) > 0) {
                        tice(i, j, k) = tice0;
                    } else {
                        tice(i, j, k) = 0.;
                    }
                }

            }
        }
    }
    advection.setData(hiceName, hice);
    advection.setData(ciceName, cice);
    advection.setData(uName, u);
    advection.setData(vName, v);

    double timeLimit = 1e5;
    TimePoint origin("2024-08-01T00:00:00Z");
    // Iterate: run the DynamicsKernel update, then advect the snow
    double outPeriod = 1000;
    double outCount = outPeriod;

    // A map of the results at the target times
    std::map<double, ModelArray> hiceTestData;
    std::map<double, ModelArray> ticeTestData;

    std::vector<double> testTimes;

    for (double t = 0; t <= timeLimit; t += deltaT) {
        if (outCount >= outPeriod) {
//            writeDataPGM(-tice * 100, std::to_string(std::lround(t)) + ".pgm", 1);
//            writeDataPGM(hsnow * 100, std::to_string(std::lround(t)) + ".pgm");
// Store the data as a standard sinusoid: subtract the offset (1.5) and normalize the range (*2)
            hiceTestData[t] = 2 * (hsnow - 1.5);
            ticeTestData[t] = 2 * (tice + 1.5);
            testTimes.push_back(t);
            outCount = 0;
        }
        TimestepTime tst = { origin + Duration(std::to_string(std::lround(t))), Duration(deltaT) };
        advection.update(tst);
        advection.advectField(hsnow, hsnowName);
        advection.advectField(tice, ticeName);
        outCount += deltaT;
    }

    // Test the collected data
    // Compare the initial data with a quarter rotation later (should sum to ~0)
    size_t pTest = x0 / 2;
    // Sample 8 points around the circle to do a basic Fourier transform
    const size_t nSamples = 8;
    std::array<size_t, nSamples> xTest = { pTest, pTest, 0, -pTest, -pTest, -pTest, 0, pTest };
    std::array<size_t, nSamples> yTest = { 0, pTest, pTest, pTest, 0, -pTest, -pTest, -pTest };

    std::vector<double> hSnowExpectedCosines = { 4 * hsnowA, -4 * hsnowA, 4 * hsnowA };
    std::vector<double> hSnowExpectedSines = { 0, 0, 0 };

    for (double testTime : testTimes) {



    }

    double hiceSumSelf = 0.;
    double hiceSumQuarter = 0.;
    double hiceSumHalf = 0.;
    size_t hiceCount = 0;
    double ticeSumSelf = 0.;
    double ticeSumQuarter = 0.;
    double ticeSumHalf = 0.;
    size_t ticeCount = 0;
    size_t testRow = ny/3;
    size_t testLevel = 1;
    REQUIRE(testTimes.size() >= 3);
    for (size_t i = 0; i < nx; ++i) {
        // -3 is the missing data value, due to the arithmetic manipulations above
        if (hiceTestData.at(testTimes[0])(i, testRow) != -3) {
            ++hiceCount;
            hiceSumSelf+= hiceTestData.at(testTimes[0])(i, testRow) + hiceTestData.at(testTimes[0])(i, testRow);
            hiceSumQuarter += hiceTestData.at(testTimes[0])(i, testRow) + hiceTestData.at(testTimes[1])(i, testRow);
            hiceSumHalf+= hiceTestData.at(testTimes[0])(i, testRow) + hiceTestData.at(testTimes[2])(i, testRow);
        }
        if (ticeTestData.at(testTimes[0])(i, testRow, testLevel) != +3) {
            ++ticeCount;
            ticeSumSelf+= ticeTestData.at(testTimes[0])(i, testRow, testLevel) + ticeTestData.at(testTimes[0])(i, testRow, testLevel);
            ticeSumQuarter += ticeTestData.at(testTimes[0])(i, testRow, testLevel) + ticeTestData.at(testTimes[1])(i, testRow, testLevel);
            ticeSumHalf+= ticeTestData.at(testTimes[0])(i, testRow, testLevel) + ticeTestData.at(testTimes[2])(i, testRow, testLevel);
        }
    }
    double eps = 1.5e-2;
    // Test the mean difference when the pattern is out of phase (quarter) and in phase (half)
    // Normalize by the sum with self and the number of data points
    // The ratio of the self sum and the sum with the in-phase pattern should be close to 1
    REQUIRE(std::fabs(hiceSumHalf / hiceSumSelf - 1) / hiceCount < eps);
    // The sum of the out-of phase pattern should be close to zero
    REQUIRE(std::fabs(hiceSumQuarter / hiceSumSelf) / hiceCount < eps);

    REQUIRE(std::fabs(ticeSumHalf / ticeSumSelf - 1) / ticeCount < eps);
    // The sum of the out-of phase pattern should be close to zero
    REQUIRE(std::fabs(ticeSumQuarter / ticeSumSelf) / ticeCount < eps);
}

}
