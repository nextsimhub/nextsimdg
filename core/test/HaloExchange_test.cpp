/*!
 * @file ModelMetadata_test.cpp
 *
 * @date 28 Jan 2025
 * @author Tom Meltzer <tdm39@cam.ac.uk>
 */

#include <doctest/extensions/doctest_mpi.h>

#include "ModelMetadata.hpp"
#include "ParaGridIO.hpp"
#include "Slice.hpp"
#include "include/DGModelArray.hpp"
#include "include/Halo.hpp"
#include "include/ModelArraySlice.hpp"
#include "include/ParametricMesh.hpp"
#include "include/dgVector.hpp"

namespace Nextsim {

const std::string testFilesDir = TEST_FILES_DIR;
const std::string file_cb = testFilesDir + "/partition_metadata_3_cb.nc";
const std::string file_pb = testFilesDir + "/partition_metadata_3_pb.nc";

static const int DG = 3;
static const bool debug = false;

// aliases
using Indexer::indexer;
using Slice = ArraySlicer::Slice;
using SliceIter = ArraySlicer::SliceIter;
using MultiDim = ArraySlicer::SliceIter::MultiDim;

// print out dgvec so that origin (0,0) is bottom left corner
static void printDGVec(DGVector<DG>& dgvec, MultiDim& dims)
{
    if (debug) {
        for (size_t j = 0; j < dims[1]; ++j) {
            for (size_t i = 0; i < dims[0]; ++i) {
                size_t pos = indexer(dims, { i, dims[1] - j - 1 });
                printf("%6.0f ", dgvec(pos, 0));
            }
            printf("\n");
        }
        printf("--------------------------------------------------\n");
    }
};

ParametricMesh initializeSmesh(size_t localNx, size_t localNy)
{
    auto CoordinateSystem = Nextsim::CARTESIAN;
    ParametricMesh smesh(CoordinateSystem);
    smesh.nx = localNx + 2;
    smesh.ny = localNy + 2;
    smesh.nnodes = smesh.nx * smesh.ny;
    smesh.nelements = smesh.nnodes;
    smesh.vertices.resize(smesh.nelements, Eigen::NoChange); // note these are uninitialised but
    // unused in this test
    return smesh;
}

void initializeHField(ModelArray& hfield, size_t localNx, size_t localNy, size_t offsetX,
    size_t offsetY, int test_rank)
{
    // initialize with mock data
    for (size_t j = 0; j < localNy; ++j) {
        for (size_t i = 0; i < localNx; ++i) {
            hfield(i, j) = (test_rank + 1) * 100 + (i + offsetX) * 10 + (j + offsetY);
        }
    }
}

void setDGVecValue(DGVector<DG>& dgvec, MultiDim meshDims, double value)
{
    for (size_t k = 0; k < DG; ++k) {
        for (size_t j = 0; j < meshDims[1]; ++j) {
            for (size_t i = 0; i < meshDims[0]; ++i) {
                size_t pos = indexer(meshDims, std::vector<size_t>({ i, j }));
                dgvec(pos, k) = value;
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

    ModelMetadata metadata(file_cb, test_comm);

    const size_t nx = metadata.globalExtentX;
    const size_t ny = metadata.globalExtentY;
    const size_t localNx = metadata.localExtentX;
    const size_t localNy = metadata.localExtentY;
    const size_t offsetX = metadata.localCornerX;
    const size_t offsetY = metadata.localCornerY;

    ModelArray::setDimension(ModelArray::Dimension::X, nx, localNx, offsetX);
    ModelArray::setDimension(ModelArray::Dimension::Y, ny, localNy, offsetY);

    // create example 2D field on each process
    auto testfield = ModelArray::HField();
    testfield.resize();
    initializeHField(testfield, localNx, localNy, offsetX, offsetY, test_rank);

    // create example 2D field for dynamics code which needs +1 halo on all edges, leading to +2 in
    // x and y
    // Create the ParametricMesh object to initialize dgvec
    ParametricMesh smesh = initializeSmesh(localNx, localNy);

    // create dgvec to simulate nextsim
    DGVector<DG> dgvec(smesh);
    MultiDim meshDims = { smesh.nx, smesh.ny };

    // initialize dgvec to -1
    setDGVecValue(dgvec, meshDims, -1.0);

    // we want to copy testfield to the central block
    Slice centreBlock({ { 1, -1 }, { 1, -1 } });
    SliceIter eigCentre(centreBlock, meshDims);

    // all spatial dims are collapsed into the first dimension (col(0))
    // The other dimensions are for higher-order DG components
    auto dgvecSpatial = dgvec.col(0);

    // copy from testfield into the dgvec
    testfield[Slice({ {}, {} })].copyToSlicedBuffer(dgvecSpatial, eigCentre);

    // print to check data copied successfully
    printDGVec(dgvec, meshDims);

    // create Halo object to facilitate exchange
    Halo halo(localNx, localNy, metadata);

    // populate send buffer (each rank does this itself)
    halo.populateSendBuffer(testfield);
    // populate recv buffer (each rank uses RMA to get data from other rank's send buffers)
    halo.populateRecvBuffer();

    // update the DGVector with the halo from the recv buffer
    halo.updateDGVec(dgvec);

    // print to check halo cells copied correctly
    printDGVec(dgvec, meshDims);

    // mock data for halo exchange for each test rank (0-2)
    std::array<std::vector<double>, 3> mockDataAllProcs = {
        std::vector<double>(
            { -1, 0, 0, 0, 0, 0, 0, 0, -1, 0, 100, 110, 120, 130, 140, 150, 160, 370, 0, 101, 111,
                121, 131, 141, 151, 161, 371, 0, 102, 112, 122, 132, 142, 152, 162, 372, 0, 103,
                113, 123, 133, 143, 153, 163, 373, -1, 204, 214, 224, 234, 244, 254, 264, -1 }),
        std::vector<double>({ -1, 103, 113, 123, 133, 143, 153, 163, -1, 0, 204, 214, 224, 234, 244,
            254, 264, 374, 0, 205, 215, 225, 235, 245, 255, 265, 375, 0, 206, 216, 226, 236, 246,
            256, 266, 376, 0, 207, 217, 227, 237, 247, 257, 267, 377, 0, 208, 218, 228, 238, 248,
            258, 268, 378, -1, 0, 0, 0, 0, 0, 0, 0, -1 }),
        std::vector<double>({ -1, 0, 0, 0, -1, 160, 370, 380, 390, 0, 161, 371, 381, 391, 0, 162,
            372, 382, 392, 0, 163, 373, 383, 393, 0, 264, 374, 384, 394, 0, 265, 375, 385, 395, 0,
            266, 376, 386, 396, 0, 267, 377, 387, 397, 0, 268, 378, 388, 398, 0, -1, 0, 0, 0, -1 }),
    };

    // create mock eigen matrix for each process
    Eigen::Matrix<double, Eigen::Dynamic, 1> mockData;
    mockData.resize(dgvecSpatial.size(), 1);
    mockData = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, 1>>(
        mockDataAllProcs[test_rank].data(), mockData.size(), 1);

    // spatial part of dgvec
    auto testData = dgvec(Eigen::all, 0);

    // compare mock data with data generated in this test (testData)
    REQUIRE(testData.isApprox(mockData));
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
    const size_t localNx = metadata.localExtentX;
    const size_t localNy = metadata.localExtentY;
    const size_t offsetX = metadata.localCornerX;
    const size_t offsetY = metadata.localCornerY;

    ModelArray::setDimension(ModelArray::Dimension::X, nx, localNx, offsetX);
    ModelArray::setDimension(ModelArray::Dimension::Y, ny, localNy, offsetY);

    // create example 2D field on each process
    auto testfield = ModelArray::HField();
    testfield.resize();
    initializeHField(testfield, localNx, localNy, offsetX, offsetY, test_rank);

    // create example 2D field for dynamics code which needs +1 halo on all edges, leading to +2 in
    // x and y
    // Create the ParametricMesh object to initialize dgvec
    ParametricMesh smesh = initializeSmesh(localNx, localNy);

    // create dgvec to simulate nextsim
    DGVector<DG> dgvec(smesh);
    MultiDim meshDims = { smesh.nx, smesh.ny };

    // initialize dgvec to -1
    setDGVecValue(dgvec, meshDims, -1.0);

    // we want to copy testfield to the central block
    Slice centreBlock({ { 1, -1 }, { 1, -1 } });
    SliceIter eigCentre(centreBlock, meshDims);

    // all spatial dims are collapsed into the first dimension (col(0))
    // The other dimensions are for higher-order DG components
    auto dgvecSpatial = dgvec.col(0);

    // copy from testfield into the dgvec
    testfield[Slice({ {}, {} })].copyToSlicedBuffer(dgvecSpatial, eigCentre);

    // print to check data copied successfully
    printDGVec(dgvec, meshDims);

    Halo halo(localNx, localNy, metadata);

    halo.populateSendBuffer(testfield);
    halo.populateRecvBuffer();
    halo.updateDGVec(dgvec);

    // print to check halo cells copied correctly
    printDGVec(dgvec, meshDims);

    // mock data for halo exchange for each test rank (0-2)
    std::array<std::vector<double>, 3> mockDataAllProcs = {
        std::vector<double>({ -1, 208, 218, 228, 238, 248, 258, 268, -1, 390, 100, 110, 120, 130,
            140, 150, 160, 370, 391, 101, 111, 121, 131, 141, 151, 161, 371, 392, 102, 112, 122,
            132, 142, 152, 162, 372, 393, 103, 113, 123, 133, 143, 153, 163, 373, -1, 204, 214, 224,
            234, 244, 254, 264, -1 }),
        std::vector<double>({ -1, 103, 113, 123, 133, 143, 153, 163, -1, 394, 204, 214, 224, 234,
            244, 254, 264, 374, 395, 205, 215, 225, 235, 245, 255, 265, 375, 396, 206, 216, 226,
            236, 246, 256, 266, 376, 397, 207, 217, 227, 237, 247, 257, 267, 377, 398, 208, 218,
            228, 238, 248, 258, 268, 378, -1, 100, 110, 120, 130, 140, 150, 160, -1 }),
        std::vector<double>({ -1, 378, 388, 398, -1, 160, 370, 380, 390, 100, 161, 371, 381, 391,
            101, 162, 372, 382, 392, 102, 163, 373, 383, 393, 103, 264, 374, 384, 394, 204, 265,
            375, 385, 395, 205, 266, 376, 386, 396, 206, 267, 377, 387, 397, 207, 268, 378, 388,
            398, 208, -1, 370, 380, 390, -1 }),
    };

    // create mock eigen matrix for each process
    Eigen::Matrix<double, Eigen::Dynamic, 1> mockData;
    mockData.resize(dgvecSpatial.size(), 1);
    mockData = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, 1>>(
        mockDataAllProcs[test_rank].data(), mockData.size(), 1);

    // spatial part of dgvec
    auto testData = dgvec(Eigen::all, 0);

    std::vector<double> v2;
    v2.resize(testData.size());
    Eigen::Matrix<double, Eigen::Dynamic, 1>::Map(&v2[0], testData.size()) = testData;

    // compare mock data with data generated in this test (testData)
    REQUIRE(testData.isApprox(mockData));
}
}
