/*!
 * @file ModelMetadata_test.cpp
 *
 * @date 19 Mar 2025
 * @author Tom Meltzer <tdm39@cam.ac.uk>
 */

#include <doctest/extensions/doctest_mpi.h>

#include "ModelMetadata.hpp"
#include "include/Halo.hpp"

namespace Nextsim {

const std::string testFilesDir = TEST_FILES_DIR;
const std::string file_cb = testFilesDir + "/partition_metadata_3_cb.nc";
const std::string file_pb = testFilesDir + "/partition_metadata_3_pb.nc";

// TODO expand test for DG > 1
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

TEST_SUITE_BEGIN("Halo exchange tests");
MPI_TEST_CASE("test halo exchange on 3 proc grid", 3)
{
    // test for a 3 proc grid
    // ┌─┬─┐
    // │1│ │
    // ├─┤2│
    // │0│ │
    // └─┴─┘ (proc id)

    ModelMetadata metadata(file_cb, test_comm);

    const size_t nx = metadata.globalExtentX;
    const size_t ny = metadata.globalExtentY;
    const size_t localNx = metadata.localExtentX + 2 * Halo::haloWidth;
    const size_t localNy = metadata.localExtentY + 2 * Halo::haloWidth;
    const size_t offsetX = metadata.localCornerX;
    const size_t offsetY = metadata.localCornerY;

    ModelArray::setDimension(ModelArray::Dimension::X, nx, localNx, offsetX);
    ModelArray::setDimension(ModelArray::Dimension::Y, ny, localNy, offsetY);

    // create example 2D field on each process
    auto testData = ModelArray::HField();
    testData.resize();

    // create halo for testData model array
    Halo halo(metadata, testData);

    // create and allocate temporary Eigen array
    ModelArray::DataType innerData;
    innerData.resize(halo.getInnerSize(), testData.nComponents());
    initializeHField(innerData, localNx - 2, localNy - 2, offsetX, offsetY, test_rank);

    // populate inner block of modelarray
    halo.populateInnerBlock(innerData);
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
MPI_TEST_CASE("test halo exchange on 3 proc grid with periodic boundary conditions", 3)
{
    // test for a 3 proc grid
    // ┌─┬─┐
    // │1│ │
    // ├─┤2│
    // │0│ │
    // └─┴─┘ (proc id)

    ModelMetadata metadata(file_pb, test_comm);

    const size_t nx = metadata.globalExtentX;
    const size_t ny = metadata.globalExtentY;
    const size_t localNx = metadata.localExtentX + 2 * Halo::haloWidth;
    const size_t localNy = metadata.localExtentY + 2 * Halo::haloWidth;
    const size_t offsetX = metadata.localCornerX;
    const size_t offsetY = metadata.localCornerY;

    ModelArray::setDimension(ModelArray::Dimension::X, nx, localNx, offsetX);
    ModelArray::setDimension(ModelArray::Dimension::Y, ny, localNy, offsetY);

    // create example 2D field on each process
    auto testData = ModelArray::HField();
    testData.resize();

    // create halo for testData model array
    Halo halo(metadata, testData);

    // create and allocate temporary Eigen array
    ModelArray::DataType innerData;
    innerData.resize(halo.getInnerSize(), testData.nComponents());
    initializeHField(innerData, localNx - 2, localNy - 2, offsetX, offsetY, test_rank);

    // populate inner block of modelarray
    halo.populateInnerBlock(innerData);
    // exchange halos
    halo.exchangeHalos();

    std::array<std::vector<double>, 3> mockDataAllProcs = {
        std::vector<double>({ 0, 208, 218, 228, 238, 248, 258, 268, 0, 390, 100, 110, 120, 130, 140,
            150, 160, 370, 391, 101, 111, 121, 131, 141, 151, 161, 371, 392, 102, 112, 122, 132,
            142, 152, 162, 372, 393, 103, 113, 123, 133, 143, 153, 163, 373, 0, 204, 214, 224, 234,
            244, 254, 264, 0 }),
        std::vector<double>({ 0, 103, 113, 123, 133, 143, 153, 163, 0, 394, 204, 214, 224, 234, 244,
            254, 264, 374, 395, 205, 215, 225, 235, 245, 255, 265, 375, 396, 206, 216, 226, 236,
            246, 256, 266, 376, 397, 207, 217, 227, 237, 247, 257, 267, 377, 398, 208, 218, 228,
            238, 248, 258, 268, 378, 0, 100, 110, 120, 130, 140, 150, 160, 0 }),
        std::vector<double>({ 0, 378, 388, 398, 0, 160, 370, 380, 390, 100, 161, 371, 381, 391, 101,
            162, 372, 382, 392, 102, 163, 373, 383, 393, 103, 264, 374, 384, 394, 204, 265, 375,
            385, 395, 205, 266, 376, 386, 396, 206, 267, 377, 387, 397, 207, 268, 378, 388, 398,
            208, 0, 370, 380, 390, 0 }),
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
