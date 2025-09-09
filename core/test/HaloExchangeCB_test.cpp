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
 */

#include <doctest/extensions/doctest_mpi.h>

#include "ModelMPI.hpp"
#include "ModelMetadata.hpp"
#include "cgVector.hpp"
#include "include/DGModelArray.hpp"
#include "include/Halo.hpp"
#include "include/Interpolations.hpp"

namespace Nextsim {

const std::string testFilesDir = TEST_FILES_DIR;
const std::string file = testFilesDir + "/partition_metadata_3_cb.nc";

static const int DG = 3;
static const int CGDEGREE = 2;
static const bool debug = false;

// reference data for each process
static std::array<std::vector<double>, 3> refDataAllProcs = {
    std::vector<double>({ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 100, 110, 120, 130, 140, 150, 160, 370, 0,
        101, 111, 121, 131, 141, 151, 161, 371, 0, 102, 112, 122, 132, 142, 152, 162, 372, 0, 103,
        113, 123, 133, 143, 153, 163, 373, 0, 204, 214, 224, 234, 244, 254, 264, 0 }),
    std::vector<double>({ 0, 103, 113, 123, 133, 143, 153, 163, 0, 0, 204, 214, 224, 234, 244, 254,
        264, 374, 0, 205, 215, 225, 235, 245, 255, 265, 375, 0, 206, 216, 226, 236, 246, 256, 266,
        376, 0, 207, 217, 227, 237, 247, 257, 267, 377, 0, 208, 218, 228, 238, 248, 258, 268, 378,
        0, 0, 0, 0, 0, 0, 0, 0, 0 }),
    std::vector<double>({ 0, 0, 0, 0, 0, 160, 370, 380, 390, 0, 161, 371, 381, 391, 0, 162, 372,
        382, 392, 0, 163, 373, 383, 393, 0, 264, 374, 384, 394, 0, 265, 375, 385, 395, 0, 266, 376,
        386, 396, 0, 267, 377, 387, 397, 0, 268, 378, 388, 398, 0, 0, 0, 0, 0, 0 }),
};

void initializeTestData(ModelArray::DataType& data, size_t localNx, size_t localNy, size_t offsetX,
    size_t offsetY, size_t numComps, int test_rank)
{
    // initialize with test data
    for (size_t j = 0; j < localNy; ++j) {
        for (size_t i = 0; i < localNx; ++i) {
            for (size_t d = 0; d < numComps; ++d) {
                data(d + i * numComps + j * localNx * numComps)
                    = (d + 1) * 1000 + (test_rank + 1) * 100 + (i + offsetX) * 10 + (j + offsetY);
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
    testData.resize();
    testData = 0.;

    // create halo for testData model array
    Halo halo(testData);

    // create and allocate temporary Eigen array
    ModelArray::DataType innerData;
    innerData.resize(halo.getInnerSize(), testData.nComponents());
    initializeTestData(
        innerData, localNx - 2, localNy - 2, offsetX, offsetY, testData.nComponents(), test_rank);

    // populate inner block of modelarray
    halo.setInnerBlock(innerData, testData.getDataRef());
    // exchange halos
    halo.exchangeHalos(testData.getDataRef());

    // initialize reference data for each proc
    ModelArray::DataType refData;
    refData.resize(localNx * localNy, 1);
    refData
        = Eigen::Map<ModelArray::DataType>(refDataAllProcs[test_rank].data(), refData.size(), 1);

    for (size_t j = 0; j < localNy; j++) {
        for (size_t i = 0; i < localNx; i++) {
            if (refData(i + j * localNx) > 0) {
                REQUIRE(testData(i, j) == 1000. + refData(i + j * localNx));
            } else {
                REQUIRE(testData(i, j) == 0.);
            }
        }
    }

    VertexField coordinates = ModelArray::VertexField();
    coordinates.resize();
    coordinates = 0.;

    Halo haloVertex(coordinates);

    // Vetex coordinates
    auto localNxVertex = localNx + 1;
    auto localNyVertex = localNy + 1;
    for (size_t i = 0; i < localNxVertex - 2 * Halo::haloWidth; ++i) {
        for (size_t j = 0; j < localNyVertex - 2 * Halo::haloWidth; ++j) {
            double x = (i + offsetX) - 0.5 - float(nx) / 2;
            double y = (j + offsetY) - 0.5 - float(ny) / 2;
            coordinates.components({ i + Halo::haloWidth, j + Halo::haloWidth })[0] = x * 100.;
            coordinates.components({ i + Halo::haloWidth, j + Halo::haloWidth })[1] = y * 100.;
        }
    }

    haloVertex.exchangeHalos(coordinates.getDataRef());

    std::vector<double> refDataVertex;

    if (test_rank == 0) {
        refDataVertex = std::vector<double>({ 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
            0., 0., 0., 0., 0., 0., 0., 0., 0., -550., -500., -450., -500., -350., -500., -250.,
            -500., -150., -500., -50., -500., 50., -500., 150., -500., 250., -500., 0., 0., -550.,
            -400., -450., -400., -350., -400., -250., -400., -150., -400., -50., -400., 50., -400.,
            150., -400., 250., -400., 0., 0., -550., -300., -450., -300., -350., -300., -250.,
            -300., -150., -300., -50., -300., 50., -300., 150., -300., 250., -300., 0., 0., -550.,
            -200., -450., -200., -350., -200., -250., -200., -150., -200., -50., -200., 50., -200.,
            150., -200., 250., -200., 0., 0., -550., -100., -450., -100., -350., -100., -250.,
            -100., -150., -100., -50., -100., 50., -100., 150., -100., 250., -100., 0., 0., -550.,
            0., -450., 0., -350., 0., -250., 0., -150., 0., -50., 0., 50., 0., 150., 0., 0., 0. });
    }
    if (test_rank == 1) {
        refDataVertex = std::vector<double>({ 0., 0., -550., -200., -450., -200., -350., -200.,
            -250., -200., -150., -200., -50., -200., 50., -200., 150., -200., 0., 0., 0., 0., -550.,
            -100., -450., -100., -350., -100., -250., -100., -150., -100., -50., -100., 50., -100.,
            150., -100., 250., -100., 0., 0., -550., 0., -450., 0., -350., 0., -250., 0., -150., 0.,
            -50., 0., 50., 0., 150., 0., 250., 0., 0., 0., -550., 100., -450., 100., -350., 100.,
            -250., 100., -150., 100., -50., 100., 50., 100., 150., 100., 250., 100., 0., 0., -550.,
            200., -450., 200., -350., 200., -250., 200., -150., 200., -50., 200., 50., 200., 150.,
            200., 250., 200., 0., 0., -550., 300., -450., 300., -350., 300., -250., 300., -150.,
            300., -50., 300., 50., 300., 150., 300., 250., 300., 0., 0., -550., 400., -450., 400.,
            -350., 400., -250., 400., -150., 400., -50., 400., 50., 400., 150., 400., 250., 400.,
            0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. });
    }
    if (test_rank == 2) {
        refDataVertex = std::vector<double>({ 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 50.,
            -500., 150., -500., 250., -500., 350., -500., 450., -500., 0., 0., 50., -400., 150.,
            -400., 250., -400., 350., -400., 450., -400., 0., 0., 50., -300., 150., -300., 250.,
            -300., 350., -300., 450., -300., 0., 0., 50., -200., 150., -200., 250., -200., 350.,
            -200., 450., -200., 0., 0., 50., -100., 150., -100., 250., -100., 350., -100., 450.,
            -100., 0., 0., 50., 0., 150., 0., 250., 0., 350., 0., 450., 0., 0., 0., 50., 100., 150.,
            100., 250., 100., 350., 100., 450., 100., 0., 0., 50., 200., 150., 200., 250., 200.,
            350., 200., 450., 200., 0., 0., 50., 300., 150., 300., 250., 300., 350., 300., 450.,
            300., 0., 0., 50., 400., 150., 400., 250., 400., 350., 400., 450., 400., 0., 0., 0., 0.,
            0., 0., 0., 0., 0., 0., 0., 0., 0., 0. });
    }

    auto coordRef = coordinates.getDataRef();
    for (size_t j = 0; j < localNyVertex; j++) {
        for (size_t i = 0; i < localNxVertex; i++) {
            for (size_t coord = 0; coord < 2; coord++) {
                REQUIRE(coordRef(i + j * localNxVertex, coord)
                    == refDataVertex[coord + i * 2 + j * localNxVertex * 2]);
            }
        }
    }
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
    testData.resize();
    testData = 0.;

    // create halo for testData model array
    Halo halo(testData);

    // create and allocate temporary Eigen array
    ModelArray::DataType innerData;
    innerData.resize(halo.getInnerSize(), testData.nComponents());
    initializeTestData(
        innerData, localNx - 2, localNy - 2, offsetX, offsetY, testData.nComponents(), test_rank);

    // populate inner block of modelarray
    halo.setInnerBlock(innerData, testData.getDataRef());
    // exchange halos
    halo.exchangeHalos(testData.getDataRef());

    // initialize reference data for each proc
    ModelArray::DataType refData;
    refData.resize(localNx * localNy, 1);
    refData
        = Eigen::Map<ModelArray::DataType>(refDataAllProcs[test_rank].data(), refData.size(), 1);

    for (size_t j = 0; j < localNy; j++) {
        for (size_t i = 0; i < localNx; i++) {
            for (size_t d = 0; d < DG; ++d) {
                if (refData(i + j * localNx) > 0) {
                    REQUIRE(testData.components({ i, j })[d]
                        == (d + 1) * 1000. + refData(i + j * localNx));
                } else {
                    REQUIRE(testData.components({ i, j })[d] == 0.);
                }
            }
        }
    }
}

MPI_TEST_CASE("DGVector and CGVector", 3)
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

    // create and allocate temporary Eigen array
    ModelArray::DataType innerData;
    innerData.resize(halo.getInnerSize(), DG);
    initializeTestData(innerData, localNx - 2, localNy - 2, offsetX, offsetY, DG, test_rank);

    // populate inner block of modelarray
    halo.setInnerBlock(innerData, testData);

    // exchange halos
    halo.exchangeHalos(testData);

    // initialize reference data for each proc
    ModelArray::DataType refData;
    refData.resize(localNx * localNy, 1);
    refData
        = Eigen::Map<ModelArray::DataType>(refDataAllProcs[test_rank].data(), refData.size(), 1);

    for (size_t j = 0; j < localNy; j++) {
        for (size_t i = 0; i < localNx; i++) {
            for (size_t d = 0; d < DG; ++d) {
                if (refData(i + j * localNx) > 0) {
                    REQUIRE(
                        testData(i + j * localNx, d) == (d + 1) * 1000. + refData(i + j * localNx));
                } else {
                    REQUIRE(testData(i + j * localNx, d) == 0.);
                }
            }
        }
    }

    // initialize CGVector (after halo exchange on DGVector)
    CGVector<CGDEGREE> cgVector;
    cgVector.resize_by_mesh(smesh);
    Interpolations::DG2CG(smesh, cgVector, testData);

    // create halo for testData model array
    Halo haloCG(testData);
    // haloCG.exchangeHalos(cgVector);
}
}
