/*!
 * @file ModelMetadata_test.cpp
 *
 * @date 18 Nov 2024
 * @author Tom Meltzer <tdm39@cam.ac.uk>
 */

#include <doctest/extensions/doctest_mpi.h>

#include "Slice.hpp"
#include "include/DGModelArray.hpp"
#include "include/ModelArraySlice.hpp"
#include "mpi.h"

// aliases
using Indexer::indexer;
using Slice = ArraySlicer::Slice;
using VBounds = ArraySlicer::Slice::VBounds;
using SliceIter = ArraySlicer::SliceIter;

static const int DG = 3;
static const bool debug = true;

namespace Nextsim {

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

TEST_SUITE_BEGIN("Halo exchange 3x1");
MPI_TEST_CASE("test halo exchange", 3)
{
    const auto nx = 6;
    const auto ny = 4;
    // simulate 3 x 1 grid
    // ┌─┬─┬─┐
    // │0│1│2│
    // └─┴─┴─┘(proc id)
    size_t localNx = nx / test_nb_procs;
    size_t offsetX = test_rank * localNx;
    ModelArray::setDimension(ModelArray::Dimension::X, nx, localNx, offsetX);
    ModelArray::setDimension(ModelArray::Dimension::Y, ny, ny, 0);

    // create example 2D field of size 2x3 (x,y) on each process
    auto testfield = ModelArray::HField();
    testfield.resize();

    // initialize with mock data
    for (size_t j = 0; j < ny; ++j) {
        for (size_t i = 0; i < localNx; ++i) {
            testfield(i, j) = (test_rank + 1) * 100 + (i + offsetX) * 10 + j;
        }
    }

    // create example 2D field for dynamics code which needs +1 halo on all edges, leading to +2 in
    // x and y i.e., 4x5 on each process
    // Create the ParametricMesh object to initialize dgvec
    auto CoordinateSystem = Nextsim::CARTESIAN;
    ParametricMesh smesh(CoordinateSystem);
    smesh.nx = localNx + 2;
    smesh.ny = ny + 2;
    smesh.nnodes = smesh.nx * smesh.ny;
    smesh.nelements = smesh.nnodes;
    smesh.vertices.resize(smesh.nelements, Eigen::NoChange);

    SliceIter::MultiDim meshDims = { smesh.nx, smesh.ny };
    for (size_t i = 0; i < nx; ++i) {
        for (size_t j = 0; j < ny; ++j) {
            size_t pos = indexer(meshDims, std::vector<size_t>({ i, j }));
            smesh.vertices(pos, 0) = i;
            smesh.vertices(pos, 1) = j;
        }
    }

    // create dgvec to simulate nextsim
    DGVector<DG> dgvec(smesh);

    // initialize dgvec to -1
    for (size_t k = 0; k < DG; ++k) {
        for (size_t j = 0; j < smesh.ny; ++j) {
            for (size_t i = 0; i < smesh.nx; ++i) {
                size_t pos = indexer(meshDims, std::vector<size_t>({ i, j }));
                dgvec(pos, k) = -1.;
            }
        }
    }

    // print to check
    printDGVec(dgvec, meshDims);

    // make some slices for later use
    // clang-format off
    Slice wholeArray {{{  },{  }}};
    Slice sliceLeft  {{{ 0},{  }}};
    Slice sliceRight {{{-1},{  }}};
    // Slice sliceTop   {{{  },{-1}}}; //TODO is this unexpected behaviour (-1 doesnt work)
    Slice sliceTop   {{{  },{ny-1}}};
    Slice sliceBottom{{{  },{ 0}}};
    // clang-format on

    // clang-format off
    // Slice wholeArray2  {{{0,localNx,1},{0,ny,1}}};
    // Slice sliceLeft2   {{{0},{0,ny,1}}};
    // Slice sliceRight2  {{{localNx},{0,ny,1}}};
    // Slice sliceTop2    {{{0,localNx,1},{ny-1}}};
    // Slice sliceBottom2 {{{0,localNx,1},{0}}};
    // clang-format on

    // we want to copy testfield to the central block
    Slice centreBlock { { { 1, -1 }, { 1, -1 } } };
    SliceIter eigCentre(centreBlock, meshDims);

    // all spatial dims are collapsed into the first dimension (col(0))
    // The other dimensions are for higher-order DG components
    auto dgvecSpatial = dgvec.array().col(0);

    // copy from testfield into the dgvec
    testfield[wholeArray].copyToSlicedBuffer(dgvecSpatial, eigCentre);

    // print to check data copied successfully
    printDGVec(dgvec, meshDims);

    int ierr;
    enum Edge { BOTTOM, RIGHT, TOP, LEFT, N_EDGE };
    constexpr std::array<Edge, N_EDGE> edges = { BOTTOM, RIGHT, TOP, LEFT };
    std::map<Edge, std::vector<size_t>> neighbours, haloStarts, halo;
    std::map<Edge, std::vector<int>> recv;
    std::map<Edge, Slice> slices = { { LEFT, sliceLeft }, { RIGHT, sliceRight }, { TOP, sliceTop },
        { BOTTOM, sliceBottom } };
    std::map<Edge, size_t> edgeLengths
        = { { LEFT, ny }, { RIGHT, ny }, { TOP, localNx }, { BOTTOM, localNx } };

    neighbours.try_emplace(BOTTOM, std::vector<size_t>({ 0, 1, 2 }));
    neighbours.try_emplace(RIGHT, std::vector<size_t>({ 1, 2, 0 }));
    neighbours.try_emplace(TOP, std::vector<size_t>({ 0, 1, 2 }));
    neighbours.try_emplace(LEFT, std::vector<size_t>({ 2, 0, 1 }));

    haloStarts.try_emplace(BOTTOM, std::vector<size_t>(3, ny + localNx));
    haloStarts.try_emplace(RIGHT, std::vector<size_t>(3, 2 * localNx + ny));
    haloStarts.try_emplace(TOP, std::vector<size_t>(3, 0));
    haloStarts.try_emplace(LEFT, std::vector<size_t>(3, localNx));

    for (auto edge : edges) {
        size_t n = edgeLengths.at(edge);
        halo.try_emplace(edge, std::vector<size_t>(3, n));
        recv.try_emplace(edge, std::vector<int>(n, 0));
    }

    // send buffer the size of the perimeter of the domain
    size_t perimeterLength = 2 * localNx + 2 * ny;
    std::vector<int> send = std::vector<int>(perimeterLength, 0);

    size_t start = 0;
    for (auto edge : edges) {
        testfield[slices.at(edge)].copyToBuffer(send, start);
        start += edgeLengths.at(edge);
    }

    MPI_Win win;
    MPI_Win_create(
        &send[0], perimeterLength * sizeof(int), sizeof(int), MPI_INFO_NULL, test_comm, &win);

    MPI_Win_fence(MPI_MODE_NOPRECEDE, win); // Fence, no preceding RMA calls

    for (auto edge : edges) {
        int fromRank = neighbours.at(edge)[test_rank];
        size_t count = halo.at(edge)[test_rank];
        size_t disp = haloStarts.at(edge)[test_rank];
        MPI_Get(&recv.at(edge)[0], count, MPI_INT, fromRank, disp, count, MPI_INT, win);
    }

    MPI_Win_fence(MPI_MODE_NOSUCCEED, win); // Fence, no RMA calls in next epoch
    MPI_Win_free(&win); // free window

    for (auto edge : edges) {
        switch (edge) {
        case BOTTOM:
            for (size_t i = 0; i < localNx; ++i) {
                size_t pos = indexer(meshDims, { i + 1, 0 });
                dgvec(pos, 0) = recv.at(edge)[i];
            }
            break;
        case TOP:
            for (size_t i = 0; i < localNx; ++i) {
                size_t pos = indexer(meshDims, { i + 1, smesh.ny - 1 });
                dgvec(pos, 0) = recv.at(edge)[i];
            }
            break;
        case LEFT:
            for (size_t j = 0; j < ny; ++j) {
                size_t pos = indexer(meshDims, { 0, j + 1 });
                dgvec(pos, 0) = recv.at(edge)[j];
            }
            break;
        case RIGHT:
            for (size_t j = 0; j < ny; ++j) {
                size_t pos = indexer(meshDims, { smesh.nx - 1, j + 1 });
                dgvec(pos, 0) = recv.at(edge)[j];
            }
            break;
        default:
            break;
        }
    }

    printDGVec(dgvec, meshDims);

    // check left neighbour
    size_t y = 1;
    int neighbourRank = neighbours.at(LEFT)[test_rank];
    REQUIRE(dgvec(indexer(meshDims, { 0, y + 1 }), 0)
        == (neighbourRank + 1) * 100.0 + ((neighbourRank + 1) * localNx - 1) * 10.0 + y);
    // check right neighbour
    y = 2;
    neighbourRank = neighbours.at(RIGHT)[test_rank];
    REQUIRE(dgvec(indexer(meshDims, { smesh.nx - 1, y + 1 }), 0)
        == (neighbourRank + 1) * 100.0 + (neighbourRank * localNx) * 10.0 + y);
    // check bottom neighbour
    size_t x = 0;
    neighbourRank = neighbours.at(BOTTOM)[test_rank];
    REQUIRE(dgvec(indexer(meshDims, { x + 1, 0 }), 0)
        == (neighbourRank + 1) * 100.0 + (neighbourRank * localNx + x) * 10.0 + 1.0 * (ny - 1));
    // check top neighbour
    x = 1;
    neighbourRank = neighbours.at(TOP)[test_rank];
    REQUIRE(dgvec(indexer(meshDims, { x + 1, ny + 1 }), 0)
        == (neighbourRank + 1) * 100.0 + (neighbourRank * localNx + x) * 10.0);
}

TEST_SUITE_BEGIN("Halo exchange 2x2");
MPI_TEST_CASE("test halo exchange", 4)
{

    const auto nx = 6;
    const auto ny = 8;
    size_t localNx = nx / 2;
    size_t localNy = ny / 2;
    size_t offsetX, offsetY;
    // simulate 2 x 2 grid
    // ┌─┬─┐
    // │2│3│
    // ├─┼─┤
    // │0│1│
    // └─┴─┘ (proc id)
    if (test_rank % 2) {
        offsetX = localNx;
    } else {
        offsetX = 0;
    }
    if (test_rank >= 2) {
        offsetY = localNy;
    } else {
        offsetY = 0;
    }
    ModelArray::setDimension(ModelArray::Dimension::X, nx, localNx, offsetX);
    ModelArray::setDimension(ModelArray::Dimension::Y, ny, localNy, offsetY);

    // create example 2D field of size 2x3 (x,y) on each process
    auto testfield = ModelArray::HField();
    testfield.resize();

    // initialize with mock data
    for (size_t j = 0; j < localNy; ++j) {
        for (size_t i = 0; i < localNx; ++i) {
            testfield(i, j) = (test_rank + 1) * 100 + (i + offsetX) * 10 + (j + offsetY);
        }
    }

    // create example 2D field for dynamics code which needs +1 halo on all edges, leading to +2 in
    // x and y i.e., 4x5 on each process
    // Create the ParametricMesh object to initialize dgvec
    auto CoordinateSystem = Nextsim::CARTESIAN;
    ParametricMesh smesh(CoordinateSystem);
    smesh.nx = localNx + 2;
    smesh.ny = localNy + 2;
    smesh.nnodes = smesh.nx * smesh.ny;
    smesh.nelements = smesh.nnodes;
    smesh.vertices.resize(smesh.nelements, Eigen::NoChange); // note these are uninitialised but

    // create dgvec to simulate nextsim
    DGVector<DG> dgvec(smesh);
    SliceIter::MultiDim meshDims = { smesh.nx, smesh.ny };

    // initialize dgvec to -1
    for (size_t k = 0; k < DG; ++k) {
        for (size_t j = 0; j < smesh.ny; ++j) {
            for (size_t i = 0; i < smesh.nx; ++i) {
                size_t pos = indexer(meshDims, std::vector<size_t>({ i, j }));
                dgvec(pos, k) = -1.;
            }
        }
    }

    // print to check
    printDGVec(dgvec, meshDims);

    // make some slices for later use
    // clang-format off
    VBounds VBAll    {{  },{         }};
    VBounds VBLeft   {{ 0},{         }};
    VBounds VBRight  {{-1},{         }};
    VBounds VBTop    {{  },{localNy-1}};
    VBounds VBBottom {{  },{        0}};
    // clang-format on

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

    int ierr;
    enum Edge { BOTTOM, RIGHT, TOP, LEFT, N_EDGE };
    constexpr std::array<Edge, N_EDGE> edges = { BOTTOM, RIGHT, TOP, LEFT };
    std::map<Edge, std::vector<size_t>> neighbours, haloStarts, halo;
    std::map<Edge, std::vector<int>> recv;
    std::map<Edge, Slice> slices
        = { { LEFT, VBLeft }, { RIGHT, VBRight }, { TOP, VBTop }, { BOTTOM, VBBottom } };
    std::map<Edge, size_t> edgeLengths
        = { { LEFT, localNy }, { RIGHT, localNy }, { TOP, localNx }, { BOTTOM, localNx } };

    // simulate periodic neighbours
    neighbours.try_emplace(BOTTOM, std::vector<size_t>({ 2, 3, 0, 1 }));
    neighbours.try_emplace(RIGHT, std::vector<size_t>({ 1, 0, 3, 2 }));
    neighbours.try_emplace(TOP, std::vector<size_t>({ 2, 3, 0, 1 }));
    neighbours.try_emplace(LEFT, std::vector<size_t>({ 1, 0, 3, 2 }));

    haloStarts.try_emplace(BOTTOM, std::vector<size_t>(4, localNy + localNx));
    haloStarts.try_emplace(RIGHT, std::vector<size_t>(4, 2 * localNx + localNy));
    haloStarts.try_emplace(TOP, std::vector<size_t>(4, 0));
    haloStarts.try_emplace(LEFT, std::vector<size_t>(4, localNx));

    for (auto edge : edges) {
        size_t n = edgeLengths.at(edge);
        halo.try_emplace(edge, std::vector<size_t>(test_nb_procs, n));
        recv.try_emplace(edge, std::vector<int>(n, 0));
    }

    // send buffer the size of the perimeter of the domain
    size_t perimeterLength = 2 * localNx + 2 * localNy;
    std::vector<int> send = std::vector<int>(perimeterLength, 0);

    size_t start = 0;
    for (auto edge : edges) {
        testfield[slices.at(edge)].copyToBuffer(send, start);
        start += edgeLengths.at(edge);
    }

    MPI_Win win;
    MPI_Win_create(
        &send[0], perimeterLength * sizeof(int), sizeof(int), MPI_INFO_NULL, test_comm, &win);

    MPI_Win_fence(MPI_MODE_NOPRECEDE, win); // Fence, no preceding RMA calls

    for (auto edge : edges) {
        int fromRank = neighbours.at(edge)[test_rank];
        size_t count = halo.at(edge)[test_rank];
        size_t disp = haloStarts.at(edge)[test_rank];
        MPI_Get(&recv.at(edge)[0], count, MPI_INT, fromRank, disp, count, MPI_INT, win);
    }

    MPI_Win_fence(MPI_MODE_NOSUCCEED, win); // Fence, no RMA calls in next epoch
    MPI_Win_free(&win); // free window

    for (auto edge : edges) {
        switch (edge) {
        case BOTTOM:
            for (size_t i = 0; i < localNx; ++i) {
                size_t pos = indexer(meshDims, { i + 1, 0 });
                dgvec(pos, 0) = recv.at(edge)[i];
            }
            break;
        case TOP:
            for (size_t i = 0; i < localNx; ++i) {
                size_t pos = indexer(meshDims, { i + 1, smesh.ny - 1 });
                dgvec(pos, 0) = recv.at(edge)[i];
            }
            break;
        case LEFT:
            for (size_t j = 0; j < localNy; ++j) {
                size_t pos = indexer(meshDims, { 0, j + 1 });
                dgvec(pos, 0) = recv.at(edge)[j];
            }
            break;
        case RIGHT:
            for (size_t j = 0; j < localNy; ++j) {
                size_t pos = indexer(meshDims, { smesh.nx - 1, j + 1 });
                dgvec(pos, 0) = recv.at(edge)[j];
            }
            break;
        default:
            break;
        }
    }

    printDGVec(dgvec, meshDims);

    // check left neighbour
    size_t y = 1;
    size_t dx = (test_rank >= 2) ? 2 : 0;
    int neighbourRank = neighbours.at(LEFT)[test_rank];
    REQUIRE(dgvec(indexer(meshDims, { 0, y + 1 }), 0)
        == (neighbourRank + 1) * 100.0 + ((neighbourRank - dx + 1) * localNx - 1) * 10.0 + y
            + offsetY);
    // check right neighbour
    y = 2;
    neighbourRank = neighbours.at(RIGHT)[test_rank];
    REQUIRE(dgvec(indexer(meshDims, { smesh.nx - 1, y + 1 }), 0)
        == (neighbourRank + 1) * 100.0 + ((neighbourRank - dx) * localNx) * 10.0 + y + offsetY);

    // check bottom neighbour
    size_t x = 0;
    dx = (test_rank >= 2) ? 0 : 2;
    neighbourRank = neighbours.at(BOTTOM)[test_rank];
    REQUIRE(dgvec(indexer(meshDims, { x + 1, 0 }), 0)
        == (neighbourRank + 1) * 100.0 + ((neighbourRank - dx) * localNx + x) * 10.0
            + 1.0 * (smesh.ny - offsetY + 1));
    // check top neighbour
    x = 1;
    neighbourRank = neighbours.at(TOP)[test_rank];
    REQUIRE(dgvec(indexer(meshDims, { x + 1, localNy + 1 }), 0)
        == (neighbourRank + 1) * 100.0 + ((neighbourRank - dx) * localNx + x) * 10.0
            + (localNy - offsetY));
}
}
