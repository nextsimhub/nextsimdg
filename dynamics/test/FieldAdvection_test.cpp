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
    double hsnow0 = 1.;
    for (size_t j = 0; j < ny; ++j) {
        for (size_t i = 0; i < nx; ++i) {
            double xx = i * dx - x0;
            double yy = j * dy - y0;
            double r = std::sqrt(xx * xx + yy * yy);
            if (r < r0) {
                hice(i, j) = hice0;
                cice(i, j) = cice0;
                u(i, j) = -omega * yy;
                v(i, j) = omega * xx;
                hsnow(i, j) = hsnow0;
            } else {
                hice(i, j) = 0.;
                cice(i, j) = 0.;
                u(i, j) = 0.;
                v(i, j) = 0.;
                hsnow(i, j) = 0.;
            }
            x(i, j) = xx;
            y(i, j) = yy;
        }
    }
    advection.setData(hiceName, hice);
    advection.setData(ciceName, cice);
    advection.setData(uName, u);
    advection.setData(vName, v);

    double hsnowMultiplier = 2;
    // Modify the data to advect
    for (size_t i = nx / 2; i < nx; ++i) {
        hsnow(i, ny / 2) *= hsnowMultiplier;
    }

    writeDataPGM(hsnow * 100, "start.pgm");

    double timeLimit = 1e4;
    TimePoint origin("2024-08-01T00:00:00Z");
    // Iterate: run the DynamicsKernel update, then advect the snow
    for (double t = 0; t < timeLimit; t += deltaT) {
        TimestepTime tst = { origin + Duration(std::to_string(std::lround(t))), Duration(deltaT) };
        advection.update(tst);
        advection.advectField(hsnow, hsnowName);
        double hsnowMax = -1.;
        size_t hsnowMaxJ = -1;
        size_t hsnowMaxI = (2 * nx) / 3;
        for (size_t j = 0; j < ny; ++j) {
            if (hsnow(hsnowMaxI, j) > hsnowMax) {
                hsnowMax = hsnow(hsnowMaxI, j);
                hsnowMaxJ = j;
            }
        }
        std::cout << t << " s: hsnowMax = " << hsnowMax << " @ j = " << hsnowMaxJ << std::endl;
        writeDataPGM(hsnow * 100, std::to_string(std::lround(t)) + ".pgm");
    }

}

}
