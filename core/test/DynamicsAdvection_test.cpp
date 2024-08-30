/*!
 * @file DynamicsAdvection_test.cpp
 *
 * @date 29 Aug 2023
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "include/AdvectionDynamicsKernel.hpp"
#include "include/DynamicsParameters.hpp"
#include "include/IDynamics.hpp"

#include <tuple>

#ifndef DGCOMP
#define DGCOMP 3 // Define to prevent errors from static analysis tools
#error "Number of DG components (DGCOMP) not defined" // But throw an error anyway
#endif

namespace Nextsim {
static const std::vector<std::string> namedFields = { hiceName, ciceName, uName, vName };

// The advection-only dynamics implementation
class AdvectionDynamics : public IDynamics {
public:
    AdvectionDynamics()
        : IDynamics()
        , kernel()
    {
    }
    std::string getName() const override { return "AdvectionDynamics"; }
    void update(const TimestepTime& tst) override
    {
        std::cout << tst.start << std::endl;

        // set the updated ice thickness and concentration
        kernel.setData(hiceName, hice.data());
        kernel.setData(ciceName, cice.data());

        kernel.update(tst);
    }
    void setData(const ModelState::DataMap& ms) override
    {
        // Degrees to radians as a hex float
        static const double radians = 0x1.1df46a2529d39p-6;

        IDynamics::setData(ms);

        bool isSpherical = checkSpherical(ms);

        ModelArray coords = ms.at(coordsName);
        if (isSpherical) {
            coords *= radians;
        }
        // TODO: Some encoding of the periodic edge boundary conditions
        kernel.initialise(coords, false, ms.at(maskName));

        // Set the data in the kernel arrays.
        for (const auto& fieldName : { hiceName, ciceName, uName, vName }) {
            kernel.setData(fieldName, ms.at(fieldName));
        }
    }

protected:
    ModelArray& advectHField(ModelArray& field, const std::string& fieldName) override
    {
        return kernel.advectField(field, fieldName);
    }

private:
    AdvectionDynamicsKernel<DGCOMP> kernel;
    DynamicsParameters params;
};

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

TEST_SUITE_BEGIN("Advection");
TEST_CASE("Advect a 2D field")
{
    // Parameters of the mesh
    const size_t nx = 101;
    const size_t ny = 101;
    const size_t nz = 3;

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
    // Create the data to initialize the model
    ModelState::DataMap state = {
            {coordsName, coords},
            {xName, x},
            {yName, y},
            {maskName, mask},
            {uName, u},
            {vName, v},
            {ciceName, cice},
            {hiceName, hice},
    };

    AdvectionDynamics advection;
    advection.setData(state);

    double timeLimit = 1e5;
    TimePoint origin("2024-08-01T00:00:00Z");
    // Iterate: run the DynamicsKernel update, then advect the snow
    double outPeriod = 25000;
    double outCount = outPeriod;
    size_t nOut = 0;

    // A map of the results at the target times
    std::map<double, std::pair<double, double>> hsnowExpdData;
    std::map<double, std::pair<double, std::array<double, nz>>> ticeExpdData;
    std::vector<double> hsnowExpdCosines = { 2 * hsnowA, -2 * hsnowA, 2 * hsnowA, -2 * hsnowA, 2 * hsnowA };
    std::vector<double> hsnowExpdSines = { 0, 0, 0, 0, 0 };
    std::vector<std::array<double, 3>> ticeExpdCosines = {
            { 2 * ticeA, 2 * ticeA, 2 * ticeA },
            { -2 * ticeA, -2 * ticeA, -2 * ticeA },
            { 2 * ticeA, 2 * ticeA, 2 * ticeA },
            { -2 * ticeA, -2 * ticeA, -2 * ticeA },
            { 2 * ticeA, 2 * ticeA, 2 * ticeA },
    };
    std::vector<std::array<double, 3>> ticeExpdSines = {
            { 0, 0, 0 },
            { 0, 0, 0 },
            { 0, 0, 0 },
            { 0, 0, 0 },
            { 0, 0, 0 },
    };

    size_t pTest = nx / 4;
    double testEps = 5e-2;

    for (double t = 0; t <= timeLimit; t += deltaT) {
        if (outCount >= outPeriod) {
//            writeDataPGM(-tice * 100, "0_" + std::to_string(std::lround(t)) + ".pgm", 0);
//            writeDataPGM(-tice * 100, "1_" + std::to_string(std::lround(t)) + ".pgm", 1);
//            writeDataPGM(hsnow * 100, std::to_string(std::lround(t)) + ".pgm");
            auto [hsnowCos, hsnowSin] = pseudoFT(hsnow, nx/2, ny/2, pTest);
            auto [t0Cos, t0Sin] = pseudoFT(tice, 0, nx/2, ny/2, pTest);
            auto [t1Cos, t1Sin] = pseudoFT(tice, 1, nx/2, ny/2, pTest);
            REQUIRE(hsnowCos == doctest::Approx(hsnowExpdCosines[nOut]).epsilon(testEps));
            REQUIRE(hsnowSin == doctest::Approx(hsnowExpdSines[nOut]).epsilon(testEps));
            REQUIRE(t0Cos == doctest::Approx(ticeExpdCosines[nOut][0]).epsilon(testEps));
            REQUIRE(t0Sin == doctest::Approx(ticeExpdSines[nOut][0]).epsilon(testEps));
            REQUIRE(t1Cos == doctest::Approx(ticeExpdCosines[nOut][1]).epsilon(testEps));
            REQUIRE(t1Sin == doctest::Approx(ticeExpdSines[nOut][1]).epsilon(testEps));
//            std::cout << "cosine: data=" << t0Cos << ", expected=" << ticeExpdCosines[nOut][1] << std::endl;
//            std::cout << "sine:   data=" << t0Sin << ", expected=" << ticeExpdSines[nOut][1] << std::endl;
            outCount = 0;
            ++nOut;
        }
        TimestepTime tst = { origin + Duration(std::to_string(std::lround(t))), Duration(deltaT) };
        advection.update(tst);
        advection.advectField(hsnow, hsnowName);
        advection.advectField(tice, ticeName);
        outCount += deltaT;
    }
}
}
