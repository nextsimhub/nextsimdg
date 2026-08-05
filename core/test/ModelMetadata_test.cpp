/*!
 * @author  Tom Meltzer <tdm39@cam.ac.uk>
 * @brief   Test ModelMetadata class for partition metadata reading
 * @details
 *
 * Test ModelMetadata class. Check that it correctly reads the partition metadata from a netCDF
 * file and populates the relevant member variables.
 *
 * The test uses a partition file generated with periodic boundary conditions (--px --py), so that
 * the unified neighbour arrays contain both the (former) non-periodic and periodic neighbours,
 * merged into the single (unsuffixed) fields and sorted by neighbour rank.
 */

#include <doctest/extensions/doctest_mpi.h>
#include <iostream>

#include "ModelMPI.hpp"
#include "ModelMetadata.hpp"

const std::string testFilesDir = TEST_FILES_DIR;
const std::string partitionFilename = testFilesDir + "/halo_pb_test_partition_metadata_3.nc";

namespace Nextsim {

using Edge = ModelMetadata::Edge;

constexpr auto BOTTOM = Edge::BOTTOM;
constexpr auto RIGHT = Edge::RIGHT;
constexpr auto TOP = Edge::TOP;
constexpr auto LEFT = Edge::LEFT;

using Corner = ModelMetadata::Corner;

constexpr auto TOP_LEFT = Corner::TOP_LEFT;
constexpr auto TOP_RIGHT = Corner::TOP_RIGHT;
constexpr auto BOTTOM_RIGHT = Corner::BOTTOM_RIGHT;
constexpr auto BOTTOM_LEFT = Corner::BOTTOM_LEFT;

typedef std::vector<int> vec;

TEST_SUITE_BEGIN("ModelMetadata");
MPI_TEST_CASE("Test getPartitionMetadata for periodic boundaries", 3)
// We only need to check periodic, because it's a superset of closed boundaries.
{
    auto& modelMPI = ModelMPI::getInstance(test_comm);
    auto& meta = ModelMetadata::getInstance(partitionFilename);
    REQUIRE(modelMPI.getComm() == test_comm);

    // this metadata is specific to the periodic boundary conditions
    if (test_rank == 0) {
        // edges
        REQUIRE(meta.neighbourRanks[LEFT] == vec { 2 });
        REQUIRE(meta.neighbourExtents[LEFT] == vec { 4 });
        REQUIRE(meta.neighbourHaloSend[LEFT] == vec { 3 });
        REQUIRE(meta.neighbourHaloRecv[LEFT] == vec { 18 });
        REQUIRE(meta.neighbourRanks[RIGHT] == vec { 2 });
        REQUIRE(meta.neighbourExtents[RIGHT] == vec { 4 });
        REQUIRE(meta.neighbourHaloSend[RIGHT] == vec { 15 });
        REQUIRE(meta.neighbourHaloRecv[RIGHT] == vec { 7 });
        REQUIRE(meta.neighbourRanks[BOTTOM] == vec { 1 });
        REQUIRE(meta.neighbourExtents[BOTTOM] == vec { 7 });
        REQUIRE(meta.neighbourHaloSend[BOTTOM] == vec { 12 });
        REQUIRE(meta.neighbourHaloRecv[BOTTOM] == vec { 0 });
        REQUIRE(meta.neighbourRanks[TOP] == vec { 1 });
        REQUIRE(meta.neighbourExtents[TOP] == vec { 7 });
        REQUIRE(meta.neighbourHaloSend[TOP] == vec { 0 });
        REQUIRE(meta.neighbourHaloRecv[TOP] == vec { 11 });

        // corners
        REQUIRE(meta.cornerRanks[TOP_LEFT] == vec { 2 });
        REQUIRE(meta.cornerHaloSend[TOP_LEFT] == vec { 7 });
        REQUIRE(meta.cornerRanks[TOP_RIGHT] == vec { 2 });
        REQUIRE(meta.cornerHaloSend[TOP_RIGHT] == vec { 19 });
        REQUIRE(meta.cornerRanks[BOTTOM_RIGHT] == vec { 2 });
        REQUIRE(meta.cornerHaloSend[BOTTOM_RIGHT] == vec { 23 });
        REQUIRE(meta.cornerRanks[BOTTOM_LEFT] == vec { 2 });
        REQUIRE(meta.cornerHaloSend[BOTTOM_LEFT] == vec { 11 });
    } else if (test_rank == 1) {
        // edges
        REQUIRE(meta.neighbourRanks[LEFT] == vec { 2 });
        REQUIRE(meta.neighbourExtents[LEFT] == vec { 5 });
        REQUIRE(meta.neighbourHaloSend[LEFT] == vec { 7 });
        REQUIRE(meta.neighbourHaloRecv[LEFT] == vec { 19 });
        REQUIRE(meta.neighbourRanks[RIGHT] == vec { 2 });
        REQUIRE(meta.neighbourExtents[RIGHT] == vec { 5 });
        REQUIRE(meta.neighbourHaloSend[RIGHT] == vec { 19 });
        REQUIRE(meta.neighbourHaloRecv[RIGHT] == vec { 7 });
        REQUIRE(meta.neighbourRanks[BOTTOM] == vec { 0 });
        REQUIRE(meta.neighbourExtents[BOTTOM] == vec { 7 });
        REQUIRE(meta.neighbourHaloSend[BOTTOM] == vec { 11 });
        REQUIRE(meta.neighbourHaloRecv[BOTTOM] == vec { 0 });
        REQUIRE(meta.neighbourRanks[TOP] == vec { 0 });
        REQUIRE(meta.neighbourExtents[TOP] == vec { 7 });
        REQUIRE(meta.neighbourHaloSend[TOP] == vec { 0 });
        REQUIRE(meta.neighbourHaloRecv[TOP] == vec { 12 });

        // corners
        REQUIRE(meta.cornerRanks[TOP_LEFT] == vec { 2 });
        REQUIRE(meta.cornerHaloSend[TOP_LEFT] == vec { 3 });
        REQUIRE(meta.cornerRanks[TOP_RIGHT] == vec { 2 });
        REQUIRE(meta.cornerHaloSend[TOP_RIGHT] == vec { 15 });
        REQUIRE(meta.cornerRanks[BOTTOM_RIGHT] == vec { 2 });
        REQUIRE(meta.cornerHaloSend[BOTTOM_RIGHT] == vec { 18 });
        REQUIRE(meta.cornerRanks[BOTTOM_LEFT] == vec { 2 });
        REQUIRE(meta.cornerHaloSend[BOTTOM_LEFT] == vec { 6 });
    } else if (test_rank == 2) {
        // edges
        REQUIRE(meta.neighbourRanks[LEFT] == vec { 0, 1 });
        REQUIRE(meta.neighbourExtents[LEFT] == vec { 4, 5 });
        REQUIRE(meta.neighbourHaloSend[LEFT] == vec { 7, 7 });
        REQUIRE(meta.neighbourHaloRecv[LEFT] == vec { 15, 19 });
        REQUIRE(meta.neighbourRanks[RIGHT] == vec { 0, 1 });
        REQUIRE(meta.neighbourExtents[RIGHT] == vec { 4, 5 });
        REQUIRE(meta.neighbourHaloSend[RIGHT] == vec { 18, 19 });
        REQUIRE(meta.neighbourHaloRecv[RIGHT] == vec { 3, 7 });
        REQUIRE(meta.neighbourRanks[BOTTOM] == vec { 2 });
        REQUIRE(meta.neighbourExtents[BOTTOM] == vec { 3 });
        REQUIRE(meta.neighbourHaloSend[BOTTOM] == vec { 12 });
        REQUIRE(meta.neighbourHaloRecv[BOTTOM] == vec { 0 });
        REQUIRE(meta.neighbourRanks[TOP] == vec { 2 });
        REQUIRE(meta.neighbourExtents[TOP] == vec { 3 });
        REQUIRE(meta.neighbourHaloSend[TOP] == vec { 0 });
        REQUIRE(meta.neighbourHaloRecv[TOP] == vec { 12 });

        // corners
        REQUIRE(meta.cornerRanks[TOP_LEFT] == vec { 0 });
        REQUIRE(meta.cornerHaloSend[TOP_LEFT] == vec { 7 });
        REQUIRE(meta.cornerRanks[TOP_RIGHT] == vec { 0 });
        REQUIRE(meta.cornerHaloSend[TOP_RIGHT] == vec { 18 });
        REQUIRE(meta.cornerRanks[BOTTOM_RIGHT] == vec { 1 });
        REQUIRE(meta.cornerHaloSend[BOTTOM_RIGHT] == vec { 23 });
        REQUIRE(meta.cornerRanks[BOTTOM_LEFT] == vec { 1 });
        REQUIRE(meta.cornerHaloSend[BOTTOM_LEFT] == vec { 11 });
    }
}

}
