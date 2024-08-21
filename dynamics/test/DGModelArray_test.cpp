/*!
 * @file DGModelArray_test.cpp
 *
 * @brief Test that the functions to convert from the dynamics code DGVector
 * to and from ModelArray function correctly.
 *
 * @date Oct 6, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "include/DGModelArray.hpp"

#include "include/ParametricMesh.hpp"

#include <list>

Nextsim::COORDINATES CoordinateSystem = Nextsim::CARTESIAN;

namespace Nextsim {

TEST_SUITE_BEGIN("DGModelArray");
TEST_CASE("DGVector from ModelArray")
{
    static const int DG = 3;
    const size_t nx = 32;
    const size_t ny = 32;
    const size_t mx = 100;
    const size_t my = 100;

    ModelArray::setDimensions(ModelArray::Type::DG, { nx, ny });
    ModelArray::setNComponents(ModelArray::Type::DG, DG);

    ModelArray source(ModelArray::Type::DG);
    source.resize();
    // Fill with data
    for (size_t i = 0; i < nx; ++i) {
        for (size_t j = 0; j < ny; j++) {
            for (size_t c = 0; c < DG; ++c) {
                source.components({ i, j })[c] = c + my * (j + mx * (i));
            }
        }
    }

    CHECK(source.components({ 14, 12 })[2] == (2 + my * (12 + mx * (14))));

    // Create the ParametricMesh object
    ParametricMesh smesh(CoordinateSystem);
    smesh.nx = nx;
    smesh.ny = ny;
    smesh.nnodes = nx * ny;
    smesh.nelements = nx * ny;
    smesh.vertices.resize(smesh.nelements, Eigen::NoChange);
    for (size_t i = 0; i < nx; ++i) {
        for (size_t j = 0; j < ny; ++j) {
            smesh.vertices(i * ny + j, 0) = i;
            smesh.vertices(i * ny + j, 1) = j;
        }
    }
    DGVector<DG> dest(smesh);
    DGModelArray::ma2dg<DG>(source, dest);

    // Did it work?
    size_t targetPoint = ModelArray::indexFromLocation(ModelArray::Type::DG, { 14, 12 });
    REQUIRE(dest(targetPoint, 2) == source.components({ 14, 12 })[2]);
}

TEST_CASE("DGVector from ModelArray::Type::H")
{
    static const int DG = 3;
    const size_t nx = 32;
    const size_t ny = 32;
    const size_t mx = 100;
    const size_t my = 100;

    ModelArray::setDimensions(ModelArray::Type::H, { nx, ny });

    ModelArray source(ModelArray::Type::H);
    source.resize();
    // Fill with data
    for (size_t i = 0; i < nx; ++i) {
        for (size_t j = 0; j < ny; j++) {
            source(i, j) = j + mx * (i);
        }
    }

    CHECK(source(14, 12) == 12 + mx * (14));

    // Create the ParametricMesh object
    ParametricMesh smesh(CoordinateSystem);
    smesh.nx = nx;
    smesh.ny = ny;
    smesh.nnodes = nx * ny;
    smesh.nelements = nx * ny;
    smesh.vertices.resize(smesh.nelements, Eigen::NoChange);
    for (size_t i = 0; i < nx; ++i) {
        for (size_t j = 0; j < ny; ++j) {
            smesh.vertices(i * ny + j, 0) = i;
            smesh.vertices(i * ny + j, 1) = j;
        }
    }
    DGVector<DG> dest(smesh);
    dest.zero(); // Zero the higher components, ma2dg does not touch them.
    DGModelArray::ma2dg<DG>(source, dest);

    // Did it work?
    size_t targetPoint = ModelArray::indexFromLocation(ModelArray::Type::H, { 14, 12 });
    REQUIRE(dest(targetPoint, 0) == source(14, 12));
    REQUIRE(dest(targetPoint, 2) == 0.);
}

TEST_CASE("ModelArray from DGVector")
{
    static const int DG = 6;
    const size_t nx = 31;
    const size_t ny = 33;
    const size_t mx = 100;
    const size_t my = 100;

    // Create the ParametricMesh object
    ParametricMesh smesh(CoordinateSystem);
    smesh.nx = nx;
    smesh.ny = ny;
    smesh.nnodes = nx * ny;
    smesh.nelements = nx * ny;
    smesh.vertices.resize(smesh.nelements, Eigen::NoChange);
    for (size_t i = 0; i < nx; ++i) {
        for (size_t j = 0; j < ny; ++j) {
            smesh.vertices(i * ny + j, 0) = i;
            smesh.vertices(i * ny + j, 1) = j;
        }
    }
    DGVector<DG> source(smesh);
    // Let the H arrays do the index calculation
    ModelArray::setDimensions(ModelArray::Type::H, { nx, ny });
    for (size_t i = 0; i < nx; ++i) {
        for (size_t j = 0; j < ny; j++) {
            size_t index = ModelArray::indexFromLocation(ModelArray::Type::H, { i, j });
            for (size_t c = 0; c < DG; ++c) {
                source(index, c) = c + my * (j + mx * (i));
            }
        }
    }
    size_t targetPoint = ModelArray::indexFromLocation(ModelArray::Type::H, { 14, 12 });
    CHECK(source(targetPoint, 4) == 4 + my * (12 + mx * (14)));

    ModelArray::setDimensions(ModelArray::Type::DG, { nx, ny });
    ModelArray::setNComponents(ModelArray::Type::DG, DG);

    ModelArray dest(ModelArray::Type::DG);
    dest.resize();
    DGModelArray::dg2ma(source, dest);

    // Did it work?
    REQUIRE(dest.components({ 14, 12 })[4] == source(targetPoint, 4));
}

TEST_CASE("Test with DG = 1") // (It would be a silly case to get wrong!)
{
    static const int DG = 1;
    const size_t nx = 23;
    const size_t ny = 27;
    const size_t mx = 100;

    // Create the ParametricMesh object
    ParametricMesh smesh(CoordinateSystem);
    smesh.nx = nx;
    smesh.ny = ny;
    smesh.nnodes = nx * ny;
    smesh.nelements = nx * ny;
    smesh.vertices.resize(smesh.nelements, Eigen::NoChange);
    for (size_t i = 0; i < nx; ++i) {
        for (size_t j = 0; j < ny; ++j) {
            smesh.vertices(i * ny + j, 0) = i;
            smesh.vertices(i * ny + j, 1) = j;
        }
    }
    DGVector<DG> source(smesh);
    // Let the H arrays do the index calculation
    ModelArray::setDimensions(ModelArray::Type::H, { nx, ny });
    for (size_t i = 0; i < nx; ++i) {
        for (size_t j = 0; j < ny; j++) {
            size_t index = ModelArray::indexFromLocation(ModelArray::Type::H, { i, j });
            for (size_t c = 0; c < DG; ++c) {
                source(index, c) = j + mx * i;
            }
        }
    }

    size_t targetPoint = ModelArray::indexFromLocation(ModelArray::Type::H, { 14, 12 });
    CHECK(source(targetPoint) == 12 + mx * 14);

    ModelArray::setDimensions(ModelArray::Type::DG, { nx, ny });
    ModelArray::setNComponents(ModelArray::Type::DG, DG);

    ModelArray dest(ModelArray::Type::DG);
    dest.resize();
    DGModelArray::dg2ma(source, dest);

    // Did it work?
    REQUIRE(dest(14, 12) == source(targetPoint));
    targetPoint = ModelArray::indexFromLocation(ModelArray::Type::H, { 12, 14 });
    REQUIRE(dest(12, 14) == source(targetPoint));

    // Change the values
    dest *= 2.;
    dest += 1.;

    targetPoint = ModelArray::indexFromLocation(ModelArray::Type::H, { 11, 13 });
    // Make sure it is different before the copy
    REQUIRE(source(targetPoint) != dest(11, 13));
    DGModelArray::ma2dg<DG>(dest, source);
    // And identical again after
    REQUIRE(source(targetPoint) == dest(11, 13));
}

// Tests to and from HFields.
TEST_CASE("Test the HField/DG0 transfer") // (It would be a silly case to get wrong!)
{
    static const int DG = 3;
    const size_t nx = 19;
    const size_t ny = 21;
    const size_t mx = 100;
    const size_t my = 100;

    // Create the ParametricMesh object
    ParametricMesh smesh(CoordinateSystem);
    smesh.nx = nx;
    smesh.ny = ny;
    smesh.nnodes = nx * ny;
    smesh.nelements = nx * ny;
    smesh.vertices.resize(smesh.nelements, Eigen::NoChange);
    for (size_t i = 0; i < nx; ++i) {
        for (size_t j = 0; j < ny; ++j) {
            smesh.vertices(i * ny + j, 0) = i;
            smesh.vertices(i * ny + j, 1) = j;
        }
    }
    DGVector<DG> source(smesh);
    // Let the H arrays do the index calculation
    ModelArray::setDimensions(ModelArray::Type::H, { nx, ny });
    for (size_t i = 0; i < nx; ++i) {
        for (size_t j = 0; j < ny; j++) {
            size_t index = ModelArray::indexFromLocation(ModelArray::Type::H, { i, j });
            for (size_t c = 0; c < DG; ++c) {
                source(index, c) = j + mx * i + mx * my * c;
            }
        }
    }

    size_t targetPoint = ModelArray::indexFromLocation(ModelArray::Type::H, { 14, 12 });
    CHECK(source(targetPoint, 0) == 12 + mx * 14);

    ModelArray::setDimensions(ModelArray::Type::DG, { nx, ny });
    ModelArray::setNComponents(ModelArray::Type::DG, DG);

    HField dest(ModelArray::Type::H);
    dest.resize();
    DGModelArray::dg2hField(source, dest);

    // Did it work?
    REQUIRE(dest(14, 12) == source(targetPoint, 0));
    targetPoint = ModelArray::indexFromLocation(ModelArray::Type::H, { 12, 14 });
    REQUIRE(dest(12, 14) == source(targetPoint, 0));

    // Change the values
    dest *= 2.;
    dest += 1.;

    targetPoint = ModelArray::indexFromLocation(ModelArray::Type::H, { 11, 13 });
    // Make sure it is different before the copy
    REQUIRE(source(targetPoint, 0) != dest(11, 13));
    double retained = source(targetPoint, 1);
    REQUIRE(source(targetPoint, 0) != retained);

    DGModelArray::hField2dg<DG>(dest, source);
    // And identical again after
    REQUIRE(source(targetPoint, 0) == dest(11, 13));
    REQUIRE(source(targetPoint, 1) == retained);
    REQUIRE(source(targetPoint, 0) != retained);
}

static const double dx = 100;
static const double dy = 100;
size_t indexr(size_t i, size_t j, size_t k)
{
    return k + dy * (j + dx * i) + 1000000;
}

TEST_CASE("Test 3D arrays")
{
    static const int DG = 3;
    const size_t nx = 13;
    const size_t ny = 17;
    const size_t nz = 3;


    ModelArray::setDimension(ModelArray::Dimension::X, nx);
    ModelArray::setDimension(ModelArray::Dimension::Y, ny);
    ModelArray::setDimension(ModelArray::Dimension::Z, nz);
    ModelArray::setNComponents(ModelArray::Type::DG, DG);

    // Create the ParametricMesh object
    ParametricMesh smesh(CoordinateSystem);
    smesh.nx = nx;
    smesh.ny = ny;
    smesh.nnodes = nx * ny;
    smesh.nelements = nx * ny;
    smesh.vertices.resize(smesh.nelements, Eigen::NoChange);
    for (size_t i = 0; i < nx; ++i) {
        for (size_t j = 0; j < ny; ++j) {
            smesh.vertices(i * ny + j, 0) = i;
            smesh.vertices(i * ny + j, 1) = j;
        }
    }

    ZField zin(ModelArray::Type::Z);
    zin.resize();
    ZField zout(ModelArray::Type::Z);
    zout.resize();
    zout = 0.;

    for (size_t k = 0; k < nz; ++k) {
        for (size_t j = 0; j < ny; ++j) {
            for (size_t i = 0; i < nx; ++i) {
                zin(i, j, k) = indexr(i, j, k);
            }
        }
    }

    for (size_t k = 0; k < nz; ++k) {
        DGVector<DG> dgData;
        dgData.resize_by_mesh(smesh);
        dgData.zero();
        DGModelArray::ma2dg<DG>(zin, dgData, k);

        // Check the corners
        REQUIRE(dgData(0, 0) == indexr(0, 0, k));
        REQUIRE(dgData(nx - 1, 0) == indexr(nx - 1, 0, k));
        REQUIRE(dgData((ny - 1) * nx, 0) == indexr(0, ny - 1, k));

        // check all of dgData
        bool isDGDataCorrect = true;
        for (size_t idx = 0; idx < nx * ny; ++idx) {
            auto loc = ModelArray::locationFromIndex(ModelArray::Type::H, idx);
            if (dgData(idx, 0) != indexr(loc[0], loc[1], k)) {
                isDGDataCorrect = false;
                std::string locString = "(" + std::to_string(loc[0]) + "," + std::to_string(loc[1]) + ")";
                REQUIRE_MESSAGE(dgData(idx, 0) != indexr(loc[0], loc[1], k), "At ", locString);
            }
        }
        if (isDGDataCorrect) MESSAGE("k=", k, ": DGData seems correct, ma2dg(ma, dg, k) succeeded");

        DGModelArray::dg2ma<DG>(dgData, zout, k);
        // Check that all higher z-levels are still zero
        double absSum = 0.;
        std::list<size_t> badElements;
        for (size_t l = k + 1; l < nz; ++l) {
            for (size_t i = 0; i < nx * ny; ++i) {
                size_t idx = i + zout.indexFromLocation({ 0UL, 0UL, l });
                double absValue = std::fabs(zout[idx]);
                if (absValue > 0) badElements.push_back(idx);
                absSum += absValue;
            }
            if (badElements.size() > 0) {
                for (auto badIdx : badElements) {
                    auto badLoc = zout.locationFromIndex(badIdx);
                    std::cout << "(";
                    for (auto badLocDim : badLoc) {
                        std::cout << badLocDim << ",";
                    }
                    std::cout << ") = " << zout(badLoc);
                    bool foundCulprit = false;
                    for (size_t idg = 0; idg < nx * ny; ++idg) {
                        if (dgData(idg, 0) == zout(badLoc)) {
                            foundCulprit = true;
                            std::cout << " = dgData(";
                            auto dgBadLoc = zout.locationFromIndex(badIdx);
                            for (auto badLocDim : dgBadLoc) {
                                std::cout << badLocDim << ",";
                            }
                            std::cout << ")" << std::endl;
                        }
                    }
                    if (!foundCulprit) {
                        std::cout << " which has no match found in dgData" << std::endl;
                    }
                }
                std::cout << std::endl;
            }
            REQUIRE_MESSAGE(absSum == 0., "k=", std::to_string(l));
        }

        REQUIRE(zout(0UL, 0UL, k) == zin(0UL, 0UL, k));
        REQUIRE(zout(nx-1, 0UL, k) == zin(nx-1, 0UL, k));
        REQUIRE(zout(0UL, ny-1, k) == zin(0UL, ny-1, k));
    }
}

TEST_SUITE_END();

}
