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

    double dx = 1000; //m
    double dy = 1000; //m
    double r0 = 40000; //m
    double deltaT = 25; //s

    // Create the grid coordinates
    ModelArray::setDimension(ModelArray::Dimension::X, nx);
    ModelArray::setDimension(ModelArray::Dimension::Y, ny);
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
    
    // A map of the results at the target times
    std::map<double, ModelArray> testData;

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
            if (r < r0) {
                hice(i, j) = hice0;
                cice(i, j) = cice0;
                u(i, j) = -omega * yij;
                v(i, j) = omega * xij;
                hsnow(i, j) = hsnow0;
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
    double outPeriod = 25000;
    double outCount = outPeriod;
    std::vector<double> testTimes;
    for (double t = 0; t <= timeLimit; t += deltaT) {
        if (outCount >= outPeriod) {
            //            writeDataPGM(tice * 100 + 255, std::to_string(std::lround(t)) + ".pgm");
//            writeDataPGM(hsnow * 100, std::to_string(std::lround(t)) + ".pgm");
// Store the data as a standard sinusoid: subtract the offset (1.5) and normalize the range (*2)
            testData[t] = 2 * (hsnow - 1.5);
            testTimes.push_back(t);
            outCount = 0;
        }
        TimestepTime tst = { origin + Duration(std::to_string(std::lround(t))), Duration(deltaT) };
        advection.update(tst);
        advection.advectField(hsnow, hsnowName);
        outCount += deltaT;
    }

    // Test the collected data
    // Compare the initial data with a quarter rotation later (should sum to ~0)
    double sumSelf = 0.;
    double sumQuarter = 0.;
    double sumHalf = 0.;
    size_t testRow = ny/3;
    size_t ptCount = 0;
    REQUIRE(testTimes.size() >= 3);
    for (size_t i = 0; i < nx; ++i) {
        // -3 is the missing data value, due to the arithmetic manipulations above
        if (testData.at(testTimes[0])(i, testRow) != -3) {
            ++ptCount;
            sumSelf+= testData.at(testTimes[0])(i, testRow) + testData.at(testTimes[0])(i, testRow);
            sumQuarter += testData.at(testTimes[0])(i, testRow) + testData.at(testTimes[1])(i, testRow);
            sumHalf+= testData.at(testTimes[0])(i, testRow) + testData.at(testTimes[2])(i, testRow);
        }
    }
    double eps = 1e-2;
    // Test the mean difference when the pattern is out of phase (quarter) and in phase (half)
    // Normalize by the sum with self and the number of data points
    // The ratio of the self sum and the sum with the in-phase pattern should be close to 1
    REQUIRE(std::fabs(sumHalf / sumSelf - 1) / ptCount < eps);
    // The sum of the out-of phase pattern should be close to zero
    REQUIRE(std::fabs(sumQuarter / sumSelf) / ptCount < eps);
}

}
