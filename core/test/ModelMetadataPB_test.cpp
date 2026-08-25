/*!
 * @author  Tom Meltzer <tdm39@cam.ac.uk>
 * @brief   Test ModelMetadata class for Periodic Boundaries (PB)
 * @details
 *
 * Test ModelMetadata class for Periodic Boundaries (PB). Check that it correctly reads the
 * partition metadata from a netCDF file and populates the relevant member variables.
 */

#include <doctest/extensions/doctest_mpi.h>
#include <iostream>

#include "ModelMPI.hpp"
#include "ModelMetadata.hpp"

const std::string testFilesDir = TEST_FILES_DIR;
const std::string partitionFilenamePB = testFilesDir + "/halo_pb_test_partition_metadata_3.nc";

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

// these tests are the same for closed boundary conditions (BC) and peridic BC
static void testNonPeriodicBC(int test_rank)
{
    auto& meta = ModelMetadata::getInstance();
    if (test_rank == 0) {
        // edges
        REQUIRE(meta.neighbourRanks[LEFT].size() == 0);
        REQUIRE(meta.neighbourRanks[RIGHT] == vec { 2 });
        REQUIRE(meta.neighbourExtents[RIGHT] == vec { 4 });
        REQUIRE(meta.neighbourHaloSend[RIGHT] == vec { 15 });
        REQUIRE(meta.neighbourHaloRecv[RIGHT] == vec { 7 });
        REQUIRE(meta.neighbourRanks[BOTTOM].size() == 0);
        REQUIRE(meta.neighbourRanks[TOP] == vec { 1 });
        REQUIRE(meta.neighbourExtents[TOP] == vec { 7 });
        REQUIRE(meta.neighbourHaloSend[TOP] == vec { 0 });
        REQUIRE(meta.neighbourHaloRecv[TOP] == vec { 11 });

        // corners
        REQUIRE(meta.cornerRanks[BOTTOM_LEFT].size() == 0);
        REQUIRE(meta.cornerRanks[BOTTOM_RIGHT].size() == 0);
        REQUIRE(meta.cornerRanks[TOP_LEFT].size() == 0);
        REQUIRE(meta.cornerRanks[TOP_RIGHT] == vec { 2 });
        REQUIRE(meta.cornerHaloSend[TOP_RIGHT] == vec { 19 });
    } else if (test_rank == 1) {
        // edges
        REQUIRE(meta.neighbourRanks[LEFT].size() == 0);
        REQUIRE(meta.neighbourRanks[RIGHT] == vec { 2 });
        REQUIRE(meta.neighbourExtents[RIGHT] == vec { 5 });
        REQUIRE(meta.neighbourHaloSend[RIGHT] == vec { 19 });
        REQUIRE(meta.neighbourHaloRecv[RIGHT] == vec { 7 });
        REQUIRE(meta.neighbourRanks[BOTTOM] == vec { 0 });
        REQUIRE(meta.neighbourExtents[BOTTOM] == vec { 7 });
        REQUIRE(meta.neighbourHaloSend[BOTTOM] == vec { 11 });
        REQUIRE(meta.neighbourHaloRecv[BOTTOM] == vec { 0 });
        REQUIRE(meta.neighbourRanks[TOP].size() == 0);

        // corners
        REQUIRE(meta.cornerRanks[BOTTOM_LEFT].size() == 0);
        REQUIRE(meta.cornerRanks[BOTTOM_RIGHT] == vec { 2 });
        REQUIRE(meta.cornerHaloSend[BOTTOM_RIGHT] == vec { 18 });
        REQUIRE(meta.cornerRanks[TOP_LEFT].size() == 0);
        REQUIRE(meta.cornerRanks[TOP_RIGHT].size() == 0);
    } else if (test_rank == 2) {
        // edges
        REQUIRE(meta.neighbourRanks[LEFT] == vec { 0, 1 });
        REQUIRE(meta.neighbourExtents[LEFT] == vec { 4, 5 });
        REQUIRE(meta.neighbourHaloSend[LEFT] == vec { 7, 7 });
        REQUIRE(meta.neighbourHaloRecv[LEFT] == vec { 15, 19 });
        REQUIRE(meta.neighbourRanks[RIGHT].size() == 0);
        REQUIRE(meta.neighbourRanks[BOTTOM].size() == 0);
        REQUIRE(meta.neighbourRanks[TOP].size() == 0);

        // corners
        REQUIRE(meta.cornerRanks[BOTTOM_LEFT].size() == 0);
        REQUIRE(meta.cornerRanks[BOTTOM_RIGHT].size() == 0);
        REQUIRE(meta.cornerRanks[TOP_LEFT].size() == 0);
        REQUIRE(meta.cornerRanks[TOP_RIGHT].size() == 0);
    }
}

TEST_SUITE_BEGIN("ModelMetadata");
MPI_TEST_CASE("Test getPartitionMetadata periodic boundary", 3)
{
    auto& modelMPI = ModelMPI::getInstance(test_comm);
    auto& meta = ModelMetadata::getInstance(partitionFilenamePB);
    REQUIRE(modelMPI.getComm() == test_comm);
    // this metadata should be identical to the Closed Boundary version so we check it again
    testNonPeriodicBC(test_rank);

    // this metadata is specific to the periodic boundary conditions
    if (test_rank == 0) {
        // edges
        REQUIRE(meta.neighbourRanksPeriodic[LEFT] == vec { 2 });
        REQUIRE(meta.neighbourExtentsPeriodic[LEFT] == vec { 4 });
        REQUIRE(meta.neighbourHaloSendPeriodic[LEFT] == vec { 3 });
        REQUIRE(meta.neighbourHaloRecvPeriodic[LEFT] == vec { 18 });
        REQUIRE(meta.neighbourRanksPeriodic[RIGHT].size() == 0);
        REQUIRE(meta.neighbourRanksPeriodic[BOTTOM] == vec { 1 });
        REQUIRE(meta.neighbourExtentsPeriodic[BOTTOM] == vec { 7 });
        REQUIRE(meta.neighbourHaloSendPeriodic[BOTTOM] == vec { 12 });
        REQUIRE(meta.neighbourHaloRecvPeriodic[BOTTOM] == vec { 0 });
        REQUIRE(meta.neighbourRanksPeriodic[TOP].size() == 0);

        // corners
        REQUIRE(meta.cornerRanksPeriodic[TOP_LEFT] == vec { 2 });
        REQUIRE(meta.cornerHaloSendPeriodic[TOP_LEFT] == vec { 7 });
        REQUIRE(meta.cornerRanksPeriodic[TOP_RIGHT].size() == 0);
        REQUIRE(meta.cornerRanksPeriodic[BOTTOM_RIGHT] == vec { 2 });
        REQUIRE(meta.cornerHaloSendPeriodic[BOTTOM_RIGHT] == vec { 23 });
        REQUIRE(meta.cornerRanksPeriodic[BOTTOM_LEFT] == vec { 2 });
        REQUIRE(meta.cornerHaloSendPeriodic[BOTTOM_LEFT] == vec { 11 });
    } else if (test_rank == 1) {
        // edges
        REQUIRE(meta.neighbourRanksPeriodic[LEFT] == vec { 2 });
        REQUIRE(meta.neighbourExtentsPeriodic[LEFT] == vec { 5 });
        REQUIRE(meta.neighbourHaloSendPeriodic[LEFT] == vec { 7 });
        REQUIRE(meta.neighbourHaloRecvPeriodic[LEFT] == vec { 19 });
        REQUIRE(meta.neighbourRanksPeriodic[RIGHT].size() == 0);
        REQUIRE(meta.neighbourRanksPeriodic[BOTTOM].size() == 0);
        REQUIRE(meta.neighbourRanksPeriodic[TOP] == vec { 0 });
        REQUIRE(meta.neighbourExtentsPeriodic[TOP] == vec { 7 });
        REQUIRE(meta.neighbourHaloSendPeriodic[TOP] == vec { 0 });
        REQUIRE(meta.neighbourHaloRecvPeriodic[TOP] == vec { 12 });

        // corners
        REQUIRE(meta.cornerRanksPeriodic[TOP_LEFT] == vec { 2 });
        REQUIRE(meta.cornerHaloSendPeriodic[TOP_LEFT] == vec { 3 });
        REQUIRE(meta.cornerRanksPeriodic[TOP_RIGHT] == vec { 2 });
        REQUIRE(meta.cornerHaloSendPeriodic[TOP_RIGHT] == vec { 15 });
        REQUIRE(meta.cornerRanksPeriodic[BOTTOM_RIGHT].size() == 0);
        REQUIRE(meta.cornerRanksPeriodic[BOTTOM_LEFT] == vec { 2 });
        REQUIRE(meta.cornerHaloSendPeriodic[BOTTOM_LEFT] == vec { 6 });
    } else if (test_rank == 2) {
        // edges
        REQUIRE(meta.neighbourRanksPeriodic[LEFT].size() == 0);
        REQUIRE(meta.neighbourRanksPeriodic[RIGHT] == vec { 0, 1 });
        REQUIRE(meta.neighbourExtentsPeriodic[RIGHT] == vec { 4, 5 });
        REQUIRE(meta.neighbourHaloSendPeriodic[RIGHT] == vec { 18, 19 });
        REQUIRE(meta.neighbourHaloRecvPeriodic[RIGHT] == vec { 3, 7 });
        REQUIRE(meta.neighbourRanksPeriodic[BOTTOM] == vec { 2 });
        REQUIRE(meta.neighbourExtentsPeriodic[BOTTOM] == vec { 3 });
        REQUIRE(meta.neighbourHaloSendPeriodic[BOTTOM] == vec { 12 });
        REQUIRE(meta.neighbourHaloRecvPeriodic[BOTTOM] == vec { 0 });
        REQUIRE(meta.neighbourRanksPeriodic[TOP] == vec { 2 });
        REQUIRE(meta.neighbourExtentsPeriodic[TOP] == vec { 3 });
        REQUIRE(meta.neighbourHaloSendPeriodic[TOP] == vec { 0 });
        REQUIRE(meta.neighbourHaloRecvPeriodic[TOP] == vec { 12 });

        // corners
        REQUIRE(meta.cornerRanksPeriodic[TOP_LEFT] == vec { 0 });
        REQUIRE(meta.cornerHaloSendPeriodic[TOP_LEFT] == vec { 7 });
        REQUIRE(meta.cornerRanksPeriodic[TOP_RIGHT] == vec { 0 });
        REQUIRE(meta.cornerHaloSendPeriodic[TOP_RIGHT] == vec { 18 });
        REQUIRE(meta.cornerRanksPeriodic[BOTTOM_RIGHT] == vec { 1 });
        REQUIRE(meta.cornerHaloSendPeriodic[BOTTOM_RIGHT] == vec { 23 });
        REQUIRE(meta.cornerRanksPeriodic[BOTTOM_LEFT] == vec { 1 });
        REQUIRE(meta.cornerHaloSendPeriodic[BOTTOM_LEFT] == vec { 11 });
    }
}

}
