/*!
 * @author  Tom Meltzer <tdm39@cam.ac.uk>
 * @brief   Test halo exchange class for Closed Boundaries (CB)
 * @details
 *
 * Test halo exchange class for Closed Boundaries (CB) on following fields:
 * - HField
 * - VertexField
 * - DGField
 * - DGVector
 * - DGVectorHolder
 * - CGVector
 */

#include <doctest/extensions/doctest_mpi.h>

#include "ModelMPI.hpp"
#include "ModelMetadata.hpp"
#include "cgVector.hpp"
#include "include/DGModelArray.hpp"
#include "include/Halo.hpp"
#include "include/Interpolations.hpp"
#include "include/dgVectorHolder.hpp"

namespace Nextsim {

const std::string testFilesDir = TEST_FILES_DIR;
const std::string file = testFilesDir + "/halo_cb_test_partition_metadata_3.nc";

static const int DG = 3;
static const int CGDEGREE = 2;
static const bool debug = false;

template <typename T>
void initializeTestData(
    T& data, size_t localNx, size_t localNy, size_t offsetX, size_t offsetY, bool isCG = false)
{
    size_t nCells = 1;
    if (isCG) {
        nCells = CGDEGREE;
    }

    const auto numComps = data.cols();

    // initialize only inner block of data (to mimic behaviour of Halo::setInnerBlock())
    for (size_t j = nCells; j < localNy - nCells; ++j) {
        for (size_t i = nCells; i < localNx - nCells; ++i) {
            for (size_t d = 0; d < numComps; ++d) {
                data(d + i * numComps + j * localNx * numComps)
                    = (d + 1) * 100 + (i - nCells + offsetX) * 10 + (j - nCells + offsetY);
            }
        }
    }
}

template <typename T>
void verifyTestData(T& data, size_t localNx, size_t localNy, size_t offsetX, size_t offsetY,
    size_t nx, size_t ny, bool isCG = false)
{
    size_t nCells = 1;
    if (isCG) {
        nCells = CGDEGREE;
    }

    const auto numComps = data.cols();
    auto globalNx = nx + 2 * Halo::haloWidth * nCells - nCells;
    auto globalNy = ny + 2 * Halo::haloWidth * nCells - nCells;

    // Verify the test data
    for (size_t j = 0; j < localNy; j++) {
        for (size_t i = 0; i < localNx; i++) {
            for (size_t d = 0; d < numComps; d++) {
                size_t globalI = i + offsetX;
                size_t globalJ = j + offsetY;

                double expectedValue;
                if (globalI < nCells || globalI >= globalNx || globalJ < nCells
                    || globalJ >= globalNy) {
                    // value at boundaries should be zero in the closed-boundary (CB) case
                    expectedValue = 0.0;
                } else {
                    // Interior points - use the formula
                    expectedValue = (d + 1) * 100 + (globalI - nCells) * 10 + (globalJ - nCells);
                }

                double actualValue = data(d + i * numComps + j * localNx * numComps);

                REQUIRE(actualValue == expectedValue);
            }
        }
    }
}

TEST_SUITE_BEGIN("Halo exchange tests");
MPI_TEST_CASE("test halo exchange on 3 proc grid", 3)
{
    // test for a 3 proc grid
    // ┌─┬─┐
    // │1│ │
    // ├─┤2│
    // │0│ │
    // └─┴─┘ (proc id)

    auto& modelMPI = ModelMPI::getInstance(test_comm);
    auto& metadata = ModelMetadata::getInstance(file);

    const size_t nx = metadata.getGlobalExtentX();
    const size_t ny = metadata.getGlobalExtentY();
    const size_t localNx = metadata.getLocalExtentX() + 2 * Halo::haloWidth;
    const size_t localNy = metadata.getLocalExtentY() + 2 * Halo::haloWidth;
    const size_t offsetX = metadata.getLocalCornerX();
    const size_t offsetY = metadata.getLocalCornerY();

    ModelArray::setDimension(ModelArray::Dimension::X, nx, localNx, offsetX);
    ModelArray::setDimension(ModelArray::Dimension::Y, ny, localNy, offsetY);
    ModelArray::setDimension(ModelArray::Dimension::XVERTEX, nx + 1, localNx + 1, offsetX);
    ModelArray::setDimension(ModelArray::Dimension::YVERTEX, ny + 1, localNy + 1, offsetY);

    // create example 2D field on each process
    auto testData = ModelArray::HField();
    testData.reinitialize();
    testData = 0.;

    // create halo for testData model array
    Halo halo(testData);

    initializeTestData(testData.getDataRef(), localNx, localNy, offsetX, offsetY);

    // exchange halos
    halo.exchangeHalos(testData.getDataRef());

    // Verify the test data after halo exchange
    verifyTestData(testData.getDataRef(), localNx, localNy, offsetX, offsetY, nx, ny);

    VertexField coordinates = ModelArray::VertexField();
    coordinates.reinitialize();
    coordinates = 0.;

    Halo haloVertex(coordinates);

    initializeTestData(coordinates.getDataRef(), localNx + 1, localNy + 1, offsetX, offsetY);

    haloVertex.exchangeHalos(coordinates.getDataRef());

    verifyTestData(
        coordinates.getDataRef(), localNx + 1, localNy + 1, offsetX, offsetY, nx + 1, ny + 1);
}

MPI_TEST_CASE("DGField", 3)
{
    auto& modelMPI = ModelMPI::getInstance();
    auto& metadata = ModelMetadata::getInstance();

    const size_t nx = metadata.getGlobalExtentX();
    const size_t ny = metadata.getGlobalExtentY();
    const size_t localNx = metadata.getLocalExtentX() + 2 * Halo::haloWidth;
    const size_t localNy = metadata.getLocalExtentY() + 2 * Halo::haloWidth;
    const size_t offsetX = metadata.getLocalCornerX();
    const size_t offsetY = metadata.getLocalCornerY();

    ModelArray::setDimension(ModelArray::Dimension::X, nx, localNx, offsetX);
    ModelArray::setDimension(ModelArray::Dimension::Y, ny, localNy, offsetY);
    ModelArray::setNComponents(ModelArray::Type::DG, DG);

    // create example 2D field on each process
    auto testData = ModelArray::DGField();
    testData.reinitialize();
    testData = 0.;

    // create halo for testData model array
    Halo halo(testData);

    initializeTestData(testData.getDataRef(), localNx, localNy, offsetX, offsetY);

    // exchange halos
    halo.exchangeHalos(testData.getDataRef());

    verifyTestData(testData.getDataRef(), localNx, localNy, offsetX, offsetY, nx, ny);
}

MPI_TEST_CASE("DGVector", 3)
{
    auto& modelMPI = ModelMPI::getInstance();
    auto& metadata = ModelMetadata::getInstance();

    const size_t nx = metadata.getGlobalExtentX();
    const size_t ny = metadata.getGlobalExtentY();
    const size_t localNx = metadata.getLocalExtentX() + 2 * Halo::haloWidth;
    const size_t localNy = metadata.getLocalExtentY() + 2 * Halo::haloWidth;
    const size_t offsetX = metadata.getLocalCornerX();
    const size_t offsetY = metadata.getLocalCornerY();

    ModelArray::setDimension(ModelArray::Dimension::X, nx, localNx, offsetX);
    ModelArray::setDimension(ModelArray::Dimension::Y, ny, localNy, offsetY);
    ModelArray::setNComponents(ModelArray::Type::DG, DG);

    ParametricMesh smesh(CARTESIAN);
    smesh.nx = localNx;
    smesh.ny = localNy;
    smesh.nelements = localNx * localNy;
    smesh.nnodes = (localNx + 1) * (localNy + 1);
    smesh.vertices.resize(smesh.nnodes, Eigen::NoChange);
    for (size_t i = 0; i < localNx + 1; ++i) {
        for (size_t j = 0; j < localNy + 1; ++j) {
            smesh.vertices(i * (localNy + 1) + j, 0) = i;
            smesh.vertices(i * (localNy + 1) + j, 1) = j;
        }
    }

    // create example DGVector
    DGVector<DG> testData(smesh);
    testData.zero();

    // create halo for testData model array
    Halo halo(testData);

    initializeTestData(testData, localNx, localNy, offsetX, offsetY);

    // exchange halos
    halo.exchangeHalos(testData);

    verifyTestData(testData, localNx, localNy, offsetX, offsetY, nx, ny);
}

MPI_TEST_CASE("DGVectorHolder", 3)
{
    auto& modelMPI = ModelMPI::getInstance();
    auto& metadata = ModelMetadata::getInstance();

    const size_t nx = metadata.getGlobalExtentX();
    const size_t ny = metadata.getGlobalExtentY();
    const size_t localNx = metadata.getLocalExtentX() + 2 * Halo::haloWidth;
    const size_t localNy = metadata.getLocalExtentY() + 2 * Halo::haloWidth;
    const size_t offsetX = metadata.getLocalCornerX();
    const size_t offsetY = metadata.getLocalCornerY();

    ModelArray::setDimension(ModelArray::Dimension::X, nx, localNx, offsetX);
    ModelArray::setDimension(ModelArray::Dimension::Y, ny, localNy, offsetY);
    ModelArray::setNComponents(ModelArray::Type::DG, DG);

    ParametricMesh smesh(CARTESIAN);
    smesh.nx = localNx;
    smesh.ny = localNy;
    smesh.nelements = localNx * localNy;
    smesh.nnodes = (localNx + 1) * (localNy + 1);
    smesh.vertices.resize(smesh.nnodes, Eigen::NoChange);
    for (size_t i = 0; i < localNx + 1; ++i) {
        for (size_t j = 0; j < localNy + 1; ++j) {
            smesh.vertices(i * (localNy + 1) + j, 0) = i;
            smesh.vertices(i * (localNy + 1) + j, 1) = j;
        }
    }

    // create example DGVector
    DGVector<DG> testData(smesh);
    testData.zero();

    DGVectorHolder<DG> testHolder;
    testHolder = DGVectorHolder<DG>(testData);

    // create halo for testData model array
    Halo halo(testHolder);

    initializeTestData(static_cast<DGVector<DG>&>(testHolder), localNx, localNy, offsetX, offsetY);

    // exchange halos
    halo.exchangeHalos(static_cast<DGVector<DG>&>(testHolder));

    verifyTestData(
        static_cast<DGVector<DG>&>(testHolder), localNx, localNy, offsetX, offsetY, nx, ny);
}

MPI_TEST_CASE("CGVector", 3)
{
    auto& modelMPI = ModelMPI::getInstance();
    auto& metadata = ModelMetadata::getInstance();

    const size_t nx = metadata.getGlobalExtentX();
    const size_t ny = metadata.getGlobalExtentY();
    const size_t localNx = metadata.getLocalExtentX() + 2 * Halo::haloWidth;
    const size_t localNy = metadata.getLocalExtentY() + 2 * Halo::haloWidth;
    const size_t offsetX = metadata.getLocalCornerX();
    const size_t offsetY = metadata.getLocalCornerY();

    ModelArray::setDimension(ModelArray::Dimension::X, nx, localNx, offsetX);
    ModelArray::setDimension(ModelArray::Dimension::Y, ny, localNy, offsetY);
    ModelArray::setNComponents(ModelArray::Type::DG, DG);

    ParametricMesh smesh(CARTESIAN);
    smesh.nx = localNx;
    smesh.ny = localNy;
    smesh.nelements = localNx * localNy;
    smesh.nnodes = (localNx + 1) * (localNy + 1);
    smesh.vertices.resize(smesh.nnodes, Eigen::NoChange);
    for (size_t i = 0; i < localNx + 1; ++i) {
        for (size_t j = 0; j < localNy + 1; ++j) {
            smesh.vertices(i * (localNy + 1) + j, 0) = i;
            smesh.vertices(i * (localNy + 1) + j, 1) = j;
        }
    }

    auto CGNx = nx * CGDEGREE + 1;
    auto CGNy = ny * CGDEGREE + 1;
    auto CGLocalNx = localNx * CGDEGREE + 1;
    auto CGLocalNy = localNy * CGDEGREE + 1;
    auto CGoffsetX = offsetX * CGDEGREE;
    auto CGoffsetY = offsetY * CGDEGREE;
    bool isCG = true;

    // create CGVector
    CGVector<CGDEGREE> cgVector;
    cgVector.resize_by_mesh(smesh);
    cgVector.zero();

    // initialize data
    initializeTestData(cgVector, CGLocalNx, CGLocalNy, CGoffsetX, CGoffsetY, isCG);

    // create halo for testData model array
    Halo halo(cgVector);
    halo.exchangeHalos(cgVector);

    verifyTestData(cgVector, CGLocalNx, CGLocalNy, CGoffsetX, CGoffsetY, CGNx, CGNy, isCG);
}
}
