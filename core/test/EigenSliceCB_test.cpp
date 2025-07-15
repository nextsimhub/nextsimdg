/*!
 * @file ModelMetadata_test.cpp
 *
 * @date 15 Jul 2025
 * @author Tom Meltzer <tdm39@cam.ac.uk>
 */

#include <doctest/extensions/doctest_mpi.h>

#include "ModelMPI.hpp"
#include "ModelMetadata.hpp"
#include "include/DGModelArray.hpp"
#include "include/EigenSlice.hpp"
#include "include/Halo.hpp"
#include "include/ParametricMesh.hpp"

namespace Nextsim {

const std::string testFilesDir = TEST_FILES_DIR;
const std::string file = testFilesDir + "/partition_metadata_3_cb.nc";

static const int DG = 3;
static const bool debug = false;

void initializeHField(ModelArray::DataType& data, size_t localNx, size_t localNy, size_t offsetX,
    size_t offsetY, int test_rank)
{
    // initialize with mock data
    for (size_t j = 0; j < localNy; ++j) {
        for (size_t i = 0; i < localNx; ++i) {
            data(i + j * localNx) = (test_rank + 1) * 100 + (i + offsetX) * 10 + (j + offsetY);
        }
    }
}

void initializeDGField(ModelArray::DataType& data, size_t localNx, size_t localNy, size_t offsetX,
    size_t offsetY, int test_rank)
{
    // initialize with mock data
    for (size_t j = 0; j < localNy; ++j) {
        for (size_t i = 0; i < localNx; ++i) {
            for (size_t k = 0; k < DG; k++) {
                data(i + j * localNx, k)
                    = (test_rank + 1) * 100 + (i + offsetX) * 10 + (j + offsetY) + 1000 * (k + 1);
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
    ModelArray::setDimension(ModelArray::Dimension::DG, DG, DG, 0);

    // create example 2D field on each process
    auto testData = ModelArray::HField();
    testData.resize();

    initializeHField(testData, localNx, localNy, offsetX, offsetY, test_rank);

    using Slice = ArraySlicer::Slice;
    Slice sl { { { 1, -1 }, { 1, -1 } } };

    EigenSlice<ModelArray::DataType> maSlice(testData, sl);

    maSlice.print();

    auto testDG = ModelArray::DGField();
    testDG.resize();

    initializeDGField(testDG, localNx, localNy, offsetX, offsetY, test_rank);

    Slice sdg { { { 1, -1 }, { 1, -1 }, { 0, -1 } } };
    EigenSlice<ModelArray::DataType> dgSlice(testDG, sdg);

    dgSlice.print();

    ParametricMesh smesh(CARTESIAN);
    smesh.nx = localNx;
    smesh.ny = localNy;
    smesh.nnodes = localNx * localNy;
    smesh.nelements = localNx * localNy;
    smesh.vertices.resize(smesh.nelements, Eigen::NoChange);
    for (size_t i = 0; i < localNx; ++i) {
        for (size_t j = 0; j < localNy; ++j) {
            smesh.vertices(i * localNy + j, 0) = i;
            smesh.vertices(i * localNy + j, 1) = j;
        }
    }
    DGVector<DG> dest(smesh);
    dest.setZero();

    EigenSlice<DGVector<DG>::EigenDGVector> dgvSlice(dest, sdg, smesh);

    dgvSlice.print();

    dgSlice.copyToSlicedBuffer(dgvSlice);

    dgvSlice.print();

    // create halo for testData model array
    Halo halo(testData);

    // create and allocate temporary Eigen array
    ModelArray::DataType innerData;
    innerData.resize(halo.getInnerSize(), testData.nComponents());
    initializeHField(innerData, localNx - 2, localNy - 2, offsetX, offsetY, test_rank);

    // populate inner block of modelarray
    halo.setInnerBlock(innerData);
    // exchange halos
    halo.exchangeHalos();

    std::array<std::vector<double>, 3> mockDataAllProcs = {
        std::vector<double>({ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 100, 110, 120, 130, 140, 150, 160, 370,
            0, 101, 111, 121, 131, 141, 151, 161, 371, 0, 102, 112, 122, 132, 142, 152, 162, 372, 0,
            103, 113, 123, 133, 143, 153, 163, 373, 0, 204, 214, 224, 234, 244, 254, 264, 0 }),
        std::vector<double>({ 0, 103, 113, 123, 133, 143, 153, 163, 0, 0, 204, 214, 224, 234, 244,
            254, 264, 374, 0, 205, 215, 225, 235, 245, 255, 265, 375, 0, 206, 216, 226, 236, 246,
            256, 266, 376, 0, 207, 217, 227, 237, 247, 257, 267, 377, 0, 208, 218, 228, 238, 248,
            258, 268, 378, 0, 0, 0, 0, 0, 0, 0, 0, 0 }),
        std::vector<double>({ 0, 0, 0, 0, 0, 160, 370, 380, 390, 0, 161, 371, 381, 391, 0, 162, 372,
            382, 392, 0, 163, 373, 383, 393, 0, 264, 374, 384, 394, 0, 265, 375, 385, 395, 0, 266,
            376, 386, 396, 0, 267, 377, 387, 397, 0, 268, 378, 388, 398, 0, 0, 0, 0, 0, 0 }),
    };

    // create mock eigen matrix for each process
    ModelArray::DataType mockData;
    mockData.resize(localNx * localNy, 1);
    mockData
        = Eigen::Map<ModelArray::DataType>(mockDataAllProcs[test_rank].data(), mockData.size(), 1);

    for (size_t j = 0; j < localNy; j++) {
        for (size_t i = 0; i < localNx; i++) {
            REQUIRE(testData(i, j) == mockData(i + j * localNx));
        }
    }
}
}
