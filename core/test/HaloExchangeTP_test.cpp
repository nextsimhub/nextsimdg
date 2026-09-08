/*!
 * @author  Mikolaj Adam Kowalski <mak60@cam.ac.uk>
 * @brief   Test halo exchange class for Tripolar grid topology (TP)
 * @details
 *
 * The purpose of this test is to verify the Halo exchange for a grid with a
 * tripolar topology.
 *
 * The tripolar topology is periodic in the X direction and folds across the
 * top (north) edge via a point reflection through the middle of the top edge.
 * The bottom (south), left (west) and right (east) edges behave as closed
 * boundaries.
 *
 * We accomplish the verification by:
 *  - Assigning each value in the grid with a unique id (hash) based on its global
 *    coordinates
 *  - We take care to not assign any values in the halo region for each rank
 *  - We then perform the halo exchange
 *  - For each rank we verify that the values stored in the field match the
 *    hash
 *
 * Test halo exchange class for Tripolar topology (TP) on following fields:
 * - HField
 * - VertexField
 * - DGField
 * - DGVector
 * - CGVector
 *
 * @note This test is currently a skeleton. The topological transformation for
 *       the tripolar fold is not yet implemented in the DomainModel helper
 *       class. The test is expected to fail until this is done.
 */

#include <doctest/extensions/doctest_mpi.h>
#include <iostream>
#include <ostream>
#include <utility>

#include "ModelMPI.hpp"
#include "ModelMetadata.hpp"
#include "cgVector.hpp"
#include "include/DGModelArray.hpp"
#include "include/Halo.hpp"
#include "include/Interpolations.hpp"

namespace Nextsim {

const std::string testFilesDir = TEST_FILES_DIR;
const std::string file = testFilesDir + "/halo_tp_test_partition_metadata_3.nc";

static const int DG = 3;
static const int CGDEGREE = 2;
static const bool debug = false;

/**
 * Helper class to group all domain-dependent calculations
 *
 * It is hardcoded for a tripolar topology:
 *  - periodic in X (west <-> east wraparound)
 *  - the top (north) edge is folded via a point reflection through the middle
 *    point of the top edge: (nx / 2, ny)
 *  - the bottom (south) edge is a closed boundary
 *
 * The interpretation of the X wraparound for cell based fields (e.g. H or DG)
 * is straightforward: on the right edge, the right neighbours of the cells are
 * the cells on the left edge of the domain.
 *
 * For point based fields there is only one caveat. The points on the right edge
 * ARE the points on the left edge. There is continuity and they need to have
 * the same identity.
 *
 * The top fold is a point reflection. A point with global coordinates (i, j)
 * above the top edge (j >= ny) is mapped to the reflected point below:
 *   (nx - 1 - i, 2 * ny - 1 - j)
 * using the symmetry point (nx / 2, ny).
 *
 * @todo Implement the topological transformation for the tripolar fold in
 *       global_coordinates() and data_point_unique_id().
 */
struct DomainModel {
    int nx;
    int ny;
    int localNx;
    int localNy;
    int offsetX;
    int offsetY;
    int haloSize;
    ModelArray::Type type;

    /**
     * @brief Construct a helper representation of the calculation domain
     *
     * @param nx X-axis size of the whole domain (number of cells in x direction)
     * @param ny Same as nx but for the Y-axis
     * @param localNx X-axis size of on the current rank (number of cells in x direction)
     * @param localNy Same as localNx but for the Y-axis
     * @param offsetX X-axis offset of the current rank. Defined as the index of the bottom left
     * cell of the patch assigned to the current rank. It excludes the halo cells
     * @param offsetY Y-axis offset. See offsetX for details.
     * @param type Type of the model array)
     */
    DomainModel(
        int nx, int ny, int localNx, int localNy, int offsetX, int offsetY, ModelArray::Type type)
        : nx(nx)
        , ny(ny)
        , localNx(localNx)
        , localNy(localNy)
        , offsetX(offsetX)
        , offsetY(offsetY)
        , haloSize(1)
        , type(type)
    {
        if (type == ModelArray::Type::VERTEX) {
            //  Vertex fields have an extra entries (e.g. 3 cells in 1D have 4 vertices)
            this->nx = nx;
            this->ny = ny + 1;
            this->localNx = localNx + 1;
            this->localNy = localNy + 1;
        } else if (type == ModelArray::Type::CG) {
            //  CG fields have an extra entries
            //  Cells are subdivided by degrees in addition to the normal vertices
            this->nx = nx * CGDEGREE;
            this->ny = ny * CGDEGREE + 1;
            this->localNx = localNx * CGDEGREE + 1;
            this->localNy = localNy * CGDEGREE + 1;
            this->offsetX = offsetX * CGDEGREE;
            this->offsetY = offsetY * CGDEGREE;
            this->haloSize = CGDEGREE;
        }
    }

    /**
     * @brief Assign a unique id to each point
     *
     * We do it by enumerating the data points in the Y, X order for each component.
     */
    FloatType data_point_unique_id(int component, int i, int j) const
    {
        const auto [globalI, globalJ] = global_coordinates(i, j);
        const int domain_size = this->nx * this->ny;

        // Next power of 10 after the grid size is the component offset
        const int component_offset = [](int size) {
            int offset = 1;
            // This may result in an infinite loop if size is too large
            // But since it is used only in test we can cut some corners
            while (offset <= size) {
                offset *= 10;
            }
            return offset;
        }(domain_size);

        return component_offset * (component + 1) + globalI * this->ny + globalJ;
    }

    /**
     * @brief Map the local coorinates of the domain to the global coordinates
     *
     * Assumes the tripolar topology:
     *  - periodic in X
     *  - top edge folded via point reflection through (nx / 2, ny)
     *  - bottom edge closed
     *
     * Here we we just map cells to cells and points to point.
     * To represent function a coordinate transformation of the values is needed
     * as well. THis is, however, applied by a separate function.
     */
    std::pair<int, int> global_coordinates(int x, int y) const
    {

        // Maps a global coordinate to the new global coordinates inside the domain
        const auto periodic_reflection = [](int global, int size) {
            // The double remainder is here to deal with the negative numbers
            // For positive single remainder if enough to map the value
            // into the range [0, size)
            //
            // For the negative we need to project them into the positive range
            // by adding the size again:
            //  E.g. for size 8 and global position -1, we need to map it to 7
            //
            // Last remainder is the fixup of the `shift` for the positive inputs
            return ((global) % size + size) % size;
        };

        const auto tripolar_point_symmetry
            = [](int globalX, int globalY, int twice_x_sym, int twice_y_sym) {
                  // Since the point symmetry may contain a half, it is not representable
                  // by an integer :-/ Fortunately we only need to know 2*p to apply
                  // point transformation formula x' = 2p - x
                  // We need to take them as arguments since may be different for CG/Vertex and DG
                  // fields
                  const auto globalX_prime = twice_x_sym - globalX;
                  const auto globalY_prime = twice_y_sym - globalY;

                  return std::make_pair(globalX_prime, globalY_prime);
              };

        auto globalX = x + this->offsetX - this->haloSize;
        auto globalY = y + this->offsetY - this->haloSize;

        // We need to use different logic for the point and cell based fields
        // For cells, we apply the tripolar fold for any cell outside the domain
        // For points we need to be careful and not duplicate identity of the points
        // We select the leftmost half of the top edge to have unique identity
        // The right half are the duplicates
        if (type == ModelArray::Type::VERTEX || type == ModelArray::Type::CG) {
            if (globalY >= ny || (globalY == (ny - 1) && globalX >= nx / 2)) {
                const auto twice_symmetry_point_x = nx;
                const auto twice_symmetry_point_y = 2 * (ny - 1);
                std::tie(globalX, globalY) = tripolar_point_symmetry(
                    globalX, globalY, twice_symmetry_point_x, twice_symmetry_point_y);
            }

        } else {
            if (globalY >= ny) {
                const auto twice_symmetry_point_x = nx - 1;
                const auto twice_symmetry_point_y = 2 * ny - 1;
                std::tie(globalX, globalY) = tripolar_point_symmetry(
                    globalX, globalY, twice_symmetry_point_x, twice_symmetry_point_y);
            }
        }

        // Apply the periodicity in X
        globalX = periodic_reflection(globalX, nx);

        return { globalX, globalY };
    }
};

template <typename T> void initializeTestData(T& data, const DomainModel& domain)
{
    // nCells is the effective size of the halo
    // The interior of the domain is [nCells, localNx - nCells) x [nCells, localNy - nCells)]
    const int nCells = domain.haloSize;
    const auto numComps = data.cols();

    // initialize only inner block of data (to mimic behaviour of Halo::setInnerBlock())
    for (int j = nCells; j < domain.localNy - nCells; ++j) {
        for (int i = nCells; i < domain.localNx - nCells; ++i) {
            for (int d = 0; d < numComps; ++d) {
                data(d + i * numComps + j * domain.localNx * numComps)
                    = domain.data_point_unique_id(d, i, j);
            }
        }
    }
}

template <typename T> void verifyTestData(T& data, const DomainModel& domain)
{
    const int nCells = domain.haloSize;
    const auto numComps = static_cast<int>(data.cols());

    // Verify the test data
    for (int j = 0; j < domain.localNy; j++) {
        for (int i = 0; i < domain.localNx; i++) {
            for (int d = 0; d < numComps; d++) {
                const auto [globalI, globalJ] = domain.global_coordinates(i, j);

                const auto expectedValue = domain.data_point_unique_id(d, i, j);
                const auto actualValue = data(d + i * numComps + j * domain.localNx * numComps);
                data(d + i * numComps + j * domain.localNx * numComps) = expectedValue;

                // REQUIRE(actualValue == expectedValue);
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

    // We work on signed integers to represent extents and grid
    // coordinates
    // Size of the entire computational domain
    const int nx = metadata.getGlobalExtentX();
    const int ny = metadata.getGlobalExtentY();

    // Size of the domain of the current rank with halo
    const int localNx = metadata.getLocalExtentX() + 2 * static_cast<int>(Halo::haloWidth);
    const int localNy = metadata.getLocalExtentY() + 2 * static_cast<int>(Halo::haloWidth);

    // Position of the lower left corner of the patch (excluding halo) in the global
    // mesh
    const int offsetX = metadata.getLocalCornerX();
    const int offsetY = metadata.getLocalCornerY();

    ModelArray::setDimension(ModelArray::Dimension::X, nx, localNx, offsetX);
    ModelArray::setDimension(ModelArray::Dimension::Y, ny, localNy, offsetY);
    ModelArray::setDimension(ModelArray::Dimension::XVERTEX, nx + 1, localNx + 1, offsetX);
    ModelArray::setDimension(ModelArray::Dimension::YVERTEX, ny + 1, localNy + 1, offsetY);

    // create example 2D field on each process
    DomainModel HFieldDomain { nx, ny, localNx, localNy, offsetX, offsetY, ModelArray::Type::H };
    auto testData = ModelArray::HField();
    testData.reinitialize();
    testData = 0.;

    // create halo for testData model array
    Halo halo(testData);

    initializeTestData(testData.getDataRef(), HFieldDomain);

    // exchange halos
    halo.exchangeHalos(testData.getDataRef());

    // Verify the test data after halo exchange
    verifyTestData(testData.getDataRef(), HFieldDomain);

    DomainModel VertexFieldDomain { nx, ny, localNx, localNy, offsetX, offsetY,
        ModelArray::Type::VERTEX };
    VertexField coordinates = ModelArray::VertexField();
    coordinates.reinitialize();
    coordinates = 0.;

    Halo haloVertex(coordinates);

    initializeTestData(coordinates.getDataRef(), VertexFieldDomain);

    haloVertex.exchangeHalos(coordinates.getDataRef());

    verifyTestData(coordinates.getDataRef(), VertexFieldDomain);
}

MPI_TEST_CASE("DGField", 3)
{
    auto& modelMPI = ModelMPI::getInstance();
    auto& metadata = ModelMetadata::getInstance();

    // We work on signed integers to represent extents and grid
    // coordinates
    // Size of the entire computational domain
    const int nx = metadata.getGlobalExtentX();
    const int ny = metadata.getGlobalExtentY();

    // Size of the domain of the current rank with halo
    const int localNx = metadata.getLocalExtentX() + 2 * static_cast<int>(Halo::haloWidth);
    const int localNy = metadata.getLocalExtentY() + 2 * static_cast<int>(Halo::haloWidth);

    // Position of the lower left corner of the patch (excluding halo) in the global
    // mesh
    const int offsetX = metadata.getLocalCornerX();
    const int offsetY = metadata.getLocalCornerY();

    ModelArray::setDimension(ModelArray::Dimension::X, nx, localNx, offsetX);
    ModelArray::setDimension(ModelArray::Dimension::Y, ny, localNy, offsetY);
    ModelArray::setNComponents(ModelArray::Type::DG, DG);

    // create example 2D field on each process
    DomainModel DGFieldDomain { nx, ny, localNx, localNy, offsetX, offsetY, ModelArray::Type::DG };
    auto testData = ModelArray::DGField();
    testData.reinitialize();
    testData = 0.;

    // create halo for testData model array
    Halo halo(testData);

    initializeTestData(testData.getDataRef(), DGFieldDomain);

    // exchange halos
    halo.exchangeHalos(testData.getDataRef());

    verifyTestData(testData.getDataRef(), DGFieldDomain);
}

MPI_TEST_CASE("DGVector", 3)
{
    auto& modelMPI = ModelMPI::getInstance();
    auto& metadata = ModelMetadata::getInstance();

    // We work on signed integers to represent extents and grid
    // coordinates
    // Size of the entire computational domain
    const int nx = metadata.getGlobalExtentX();
    const int ny = metadata.getGlobalExtentY();

    // Size of the domain of the current rank with halo
    const int localNx = metadata.getLocalExtentX() + 2 * static_cast<int>(Halo::haloWidth);
    const int localNy = metadata.getLocalExtentY() + 2 * static_cast<int>(Halo::haloWidth);

    // Position of the lower left corner of the patch (excluding halo) in the global
    // mesh
    const int offsetX = metadata.getLocalCornerX();
    const int offsetY = metadata.getLocalCornerY();

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
    // The DGVector is not a ModelArray, but it is cell-centred with DG components,
    // so in terms of topology it behaves exactly like the DG field
    DomainModel DGVectorDomain { nx, ny, localNx, localNy, offsetX, offsetY, ModelArray::Type::DG };
    DGVector<DG> testData(smesh);
    testData.zero();

    // create halo for testData model array
    Halo halo(testData);

    initializeTestData(testData, DGVectorDomain);

    // exchange halos
    halo.exchangeHalos(testData);

    verifyTestData(testData, DGVectorDomain);
}

MPI_TEST_CASE("CGVector", 3)
{
    auto& modelMPI = ModelMPI::getInstance();
    auto& metadata = ModelMetadata::getInstance();

    // We work on signed integers to represent extents and grid
    // coordinates
    // Size of the entire computational domain
    const int nx = metadata.getGlobalExtentX();
    const int ny = metadata.getGlobalExtentY();

    // Size of the domain of the current rank with halo
    const int localNx = metadata.getLocalExtentX() + 2 * static_cast<int>(Halo::haloWidth);
    const int localNy = metadata.getLocalExtentY() + 2 * static_cast<int>(Halo::haloWidth);

    // Position of the lower left corner of the patch (excluding halo) in the global
    // mesh
    const int offsetX = metadata.getLocalCornerX();
    const int offsetY = metadata.getLocalCornerY();

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

    // create CGVector
    // The CGVector is not a ModelArray, but it lives on a shared-node lattice that
    // subdivides each cell CGDEGREE times, which Type::CG encodes in DomainModel
    DomainModel CGDomain { nx, ny, localNx, localNy, offsetX, offsetY, ModelArray::Type::CG };
    CGVector<CGDEGREE> cgVector;
    cgVector.resize_by_mesh(smesh);
    cgVector.zero();

    // initialize data
    initializeTestData(cgVector, CGDomain);

    // create halo and exchange
    Halo halo(cgVector);
    halo.exchangeHalos(cgVector);

    verifyTestData(cgVector, CGDomain);
}
}
