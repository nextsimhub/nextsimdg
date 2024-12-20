/*!
 * @file ModelMetadata_test.cpp
 *
 * @date 27 Jan 2025
 * @author Tom Meltzer <tdm39@cam.ac.uk>
 */

#include <bits/types/locale_t.h>
#include <cstddef>
#include <cstdlib>
#include <doctest/extensions/doctest_mpi.h>

#include <iostream>
#include <numeric>

#include "ModelMetadata.hpp"
#include "ParaGridIO.hpp"
#include "Slice.hpp"
#include "include/DGModelArray.hpp"
#include "include/ModelArraySlice.hpp"
#include "include/ParametricMesh.hpp"
#include "include/dgVector.hpp"
#include "mpi.h"

namespace Nextsim {

const std::string testFilesDir = TEST_FILES_DIR;
const std::string file_cb = testFilesDir + "/partition_metadata_3_cb.nc";
const std::string file_pb = testFilesDir + "/partition_metadata_3_pb.nc";

static const int DG = 3;
static const bool debug = false;

// aliases
using Indexer::indexer;
using Slice = ArraySlicer::Slice;
using VBounds = ArraySlicer::Slice::VBounds;
using SliceIter = ArraySlicer::SliceIter;

constexpr ModelMetadata::Edge BOTTOM = ModelMetadata::Edge::BOTTOM;
constexpr ModelMetadata::Edge RIGHT = ModelMetadata::Edge::RIGHT;
constexpr ModelMetadata::Edge TOP = ModelMetadata::Edge::TOP;
constexpr ModelMetadata::Edge LEFT = ModelMetadata::Edge::LEFT;

// make some slices for later use
// clang-format off
VBounds VBAll    {{  },{  }};
VBounds VBLeft   {{ 0},{  }};
VBounds VBRight  {{-1},{  }};
VBounds VBTop    {{  },{-1}};
VBounds VBBottom {{  },{ 0}};
// clang-format on
std::map<ModelMetadata::Edge, Slice> slices = {
    { LEFT, VBLeft },
    { RIGHT, VBRight },
    { TOP, VBTop },
    { BOTTOM, VBBottom },
};

// print out dgvec so that origin (0,0) is bottom left corner
static void printDGVec(DGVector<DG>& dgvec, SliceIter::MultiDim& dims)
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

void setDGVecValue(DGVector<DG>& dgvec, SliceIter::MultiDim meshDims, double value)
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

void updateDGVec(DGVector<DG>& dgvec, const std::vector<double>& recv,
    const SliceIter::MultiDim meshDims, std::array<size_t, 4>& edgeLengths,
    const ModelMetadata::Edge edge)
{
    SliceIter sliceIter = SliceIter(slices.at(edge), meshDims);
    size_t offset = std::accumulate(edgeLengths.begin(), edgeLengths.begin() + edge, 0);
    std::vector<size_t> edgeIndices;
    while (!sliceIter.isEnd()) {
        const size_t start = sliceIter.index();
        const size_t step = sliceIter.step(0);
        const size_t n = sliceIter.nElements(0);
        for (int i = 0; i < n; ++i) {
            auto idx = start + i * step;
            edgeIndices.push_back(idx);
        }
        sliceIter.incrementDim(1);
    }

    for (size_t i = 0; i < edgeIndices.size() - 2; ++i) {
        dgvec(edgeIndices[i + 1], 0) = recv[offset + i];
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
    SliceIter::MultiDim meshDims = { smesh.nx, smesh.ny };

    // initialize dgvec to -1
    setDGVecValue(dgvec, meshDims, -1.0);

    // we want to copy testfield to the central block
    Slice centreBlock(VBounds({ { 1, -1 }, { 1, -1 } }));
    SliceIter eigCentre(centreBlock, meshDims);

    // all spatial dims are collapsed into the first dimension (col(0))
    // The other dimensions are for higher-order DG components
    auto dgvecSpatial = dgvec.array().col(0);

    // copy from testfield into the dgvec
    testfield[Slice(VBAll)].copyToSlicedBuffer(dgvecSpatial, eigCentre);

    // print to check data copied successfully
    printDGVec(dgvec, meshDims);

    auto edges = ModelMetadata::edges;
    std::array<size_t, 4> edgeLengths = { localNx, localNy, localNx, localNy }; // BRTL

    // create a send buffer the size of the perimeter of the domain
    // each process will populate the send buffer with their perimeter cells
    const size_t perimeterLength = 2 * localNx + 2 * localNy;
    std::vector<double> send = std::vector<double>(perimeterLength, 0.0);
    std::vector<double> recv = std::vector<double>(perimeterLength, 0.0);

    for (auto edge : edges) {
        size_t offset = std::accumulate(edgeLengths.begin(), edgeLengths.begin() + edge, 0);
        testfield[slices.at(edge)].copyToBuffer(send, offset);
    }

    // create a RMA memory window which all process will be able to access
    MPI_Win win;
    MPI_Win_create(
        &send[0], perimeterLength * sizeof(double), sizeof(double), MPI_INFO_NULL, test_comm, &win);

    MPI_Win_fence(MPI_MODE_NOPRECEDE, win); // Fence, no preceding RMA calls

    // Each process can now gather the data it needs (based on mpi metadata) from the respective
    // send buffers on the other processes (or current process under some circumstances e.g.,
    // periodic boundary conditions)

    // get non-periodic neighbours and populate recv buffer
    for (auto edge : edges) {
        auto numNeighbours = metadata.neighbourRanks[edge].size();
        if (numNeighbours) {
            // get data for each neighbour
            for (size_t i = 0; i < numNeighbours; ++i) {
                int fromRank = metadata.neighbourRanks[edge][i];
                size_t count = metadata.neighbourExtents[edge][i];
                size_t disp = metadata.neighbourHaloSend[edge][i];
                size_t recvOffset = metadata.neighbourHaloRecv[edge][i];
                MPI_Get(
                    &recv[recvOffset], count, MPI_DOUBLE, fromRank, disp, count, MPI_DOUBLE, win);
            }
        }
    }

    MPI_Win_fence(MPI_MODE_NOSUCCEED, win); // Fence, no RMA calls in next epoch
    MPI_Win_free(&win); // free window

    // populate dgvec with halo data from the recv buffer
    for (auto edge : edges) {
        updateDGVec(dgvec, recv, meshDims, edgeLengths, edge);
    }

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
    SliceIter::MultiDim meshDims = { smesh.nx, smesh.ny };

    // initialize dgvec to -1
    setDGVecValue(dgvec, meshDims, -1.0);

    // we want to copy testfield to the central block
    Slice centreBlock(VBounds({ { 1, -1 }, { 1, -1 } }));
    SliceIter eigCentre(centreBlock, meshDims);

    // all spatial dims are collapsed into the first dimension (col(0))
    // The other dimensions are for higher-order DG components
    auto dgvecSpatial = dgvec.array().col(0);

    // copy from testfield into the dgvec
    testfield[Slice(VBAll)].copyToSlicedBuffer(dgvecSpatial, eigCentre);

    // print to check data copied successfully
    printDGVec(dgvec, meshDims);

    auto edges = ModelMetadata::edges;
    std::array<size_t, 4> edgeLengths = { localNx, localNy, localNx, localNy }; // BRTL

    // create a send buffer the size of the perimeter of the domain
    // each process will populate the send buffer with their perimeter cells
    const size_t perimeterLength = 2 * localNx + 2 * localNy;
    std::vector<double> send = std::vector<double>(perimeterLength, 0.0);
    std::vector<double> recv = std::vector<double>(perimeterLength, 0.0);

    for (auto edge : edges) {
        size_t offset = std::accumulate(edgeLengths.begin(), edgeLengths.begin() + edge, 0);
        testfield[slices.at(edge)].copyToBuffer(send, offset);
    }

    // create a RMA memory window which all process will be able to access
    MPI_Win win;
    MPI_Win_create(
        &send[0], perimeterLength * sizeof(double), sizeof(double), MPI_INFO_NULL, test_comm, &win);

    MPI_Win_fence(MPI_MODE_NOPRECEDE, win); // Fence, no preceding RMA calls

    // Each process can now gather the data it needs (based on mpi metadata) from the respective
    // send buffers on the other processes (or current process under some circumstances e.g.,
    // periodic boundary conditions)

    // get non-periodic neighbours and populate recv buffer
    for (auto edge : edges) {
        auto numNeighbours = metadata.neighbourRanks[edge].size();
        if (numNeighbours) {
            // get data for each neighbour
            for (size_t i = 0; i < numNeighbours; ++i) {
                int fromRank = metadata.neighbourRanks[edge][i];
                size_t count = metadata.neighbourExtents[edge][i];
                size_t disp = metadata.neighbourHaloSend[edge][i];
                size_t recvOffset = metadata.neighbourHaloRecv[edge][i];
                MPI_Get(
                    &recv[recvOffset], count, MPI_DOUBLE, fromRank, disp, count, MPI_DOUBLE, win);
            }
        }
    }

    // get periodic neighbours and populate recv buffer
    for (auto edge : edges) {
        auto numNeighbours = metadata.neighbourRanksPeriodic[edge].size();
        if (numNeighbours) {
            // get data for each neighbour
            for (size_t i = 0; i < numNeighbours; ++i) {
                int fromRank = metadata.neighbourRanksPeriodic[edge][i];
                size_t count = metadata.neighbourExtentsPeriodic[edge][i];
                size_t disp = metadata.neighbourHaloSendPeriodic[edge][i];
                size_t recvOffset = metadata.neighbourHaloRecvPeriodic[edge][i];
                MPI_Get(
                    &recv[recvOffset], count, MPI_DOUBLE, fromRank, disp, count, MPI_DOUBLE, win);
            }
        }
    }

    MPI_Win_fence(MPI_MODE_NOSUCCEED, win); // Fence, no RMA calls in next epoch
    MPI_Win_free(&win); // free window

    // populate dgvec with halo data from the recv buffer
    for (auto edge : edges) {
        updateDGVec(dgvec, recv, meshDims, edgeLengths, edge);
    }

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
