/*!
 * @file FieldAdvection_test.cpp
 *
 * @date Aug 29, 2024
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "include/AdvectionDynamicsKernel.hpp"
#include "include/ModelArray.hpp"
#include "include/Time.hpp"

#include <fstream>

namespace Nextsim {

// Analyse the data by a hand-rolled, one frequency Fourier transform.
const static size_t nSamples = 8;

std::pair<double, double> pseudoFTKernel(double data, size_t eighth) {
    std::array<double, nSamples> cosFactor = { 1, 0, -1, 0, 1, 0, -1, 0 };
    std::array<double, nSamples> sinFactor = { 0, 1, 0, -1, 0, 1, 0, -1 };

    return {cosFactor[eighth] * data, sinFactor[eighth] * data};
}

size_t ftx(size_t x0, size_t r, size_t eighth)
{
    const std::array<int, nSamples> xOffset = { +1, +1, 0, -1, -1, -1, 0, +1 };
    return x0 + r * xOffset[eighth];
}

size_t fty(size_t y0, size_t r, size_t eighth)
{
    const std::array<int, nSamples> yOffset = { 0, +1, +1, +1, 0, -1, -1, -1 };
    return y0 + r * yOffset[eighth];
}

std::pair<double, double> pseudoFT(const ModelArray& data, size_t x0, size_t y0, size_t r)
{
    // Sample 8 points around the circle to do a basic Fourier transform

    double cosTerm = 0;
    double sinTerm = 0;
    for (size_t theta = 0; theta < nSamples; ++theta) {
        auto complexTerm = pseudoFTKernel(data(ftx(x0, r, theta), fty(x0, r, theta)), theta);
        cosTerm += complexTerm.first;
        sinTerm += complexTerm.second;
    }
    return {cosTerm, sinTerm};
}

std::pair<double, double> pseudoFT(const ModelArray& data, size_t k, size_t x0, size_t y0, size_t r)
{
    // Sample 8 points around the circle to do a basic Fourier transform

    double cosTerm = 0;
    double sinTerm = 0;
    for (size_t theta = 0; theta < nSamples; ++theta) {
        auto complexTerm = pseudoFTKernel(data(ftx(x0, r, theta), fty(x0, r, theta), k), theta);
        cosTerm += complexTerm.first;
        sinTerm += complexTerm.second;
    }
    return {cosTerm, sinTerm};
}

TEST_SUITE_BEGIN("Field advection");
TEST_CASE("Advect a field")
{
    // Parameters of the mesh
    const size_t nx = 51;
    const size_t ny = 51;
    const size_t nz = 3;

    double dx = 1000; //m
    double dy = 1000; //m
    double r0 = 20 * dx; //m
    double deltaT = 100; //s

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
    size_t nOut = 0;

    // A map of the results at the target times
    std::map<double, std::pair<double, double>> hsnowExpdData;
    std::vector<double> hsnowExpdCosines = { 2 * hsnowA, -2 * hsnowA, 2 * hsnowA, -2 * hsnowA, 2 * hsnowA };
    std::vector<double> hsnowExpdSines = { 0, 0, 0, 0, 0 };

    size_t pTest = nx / 4;
    double testEps = 5e-2;

    for (double t = 0; t <= timeLimit; t += deltaT) {
        if (outCount >= outPeriod) {
            auto [hsnowCos, hsnowSin] = pseudoFT(hsnow, nx/2, ny/2, pTest);
            REQUIRE(hsnowCos == doctest::Approx(hsnowExpdCosines[nOut]).epsilon(testEps));
            REQUIRE(hsnowSin == doctest::Approx(hsnowExpdSines[nOut]).epsilon(testEps));
            outCount = 0;
            ++nOut;
        }
        outCount += deltaT;
        TimestepTime tst = { origin + Duration(std::to_string(std::lround(t))), Duration(deltaT) };
        advection.update(tst);
        advection.advectField(hsnow, hsnowName);
    }
}

}
