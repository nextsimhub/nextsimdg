/*!
 * @author  Tom Meltzer <tdm39@cam.ac.uk>
 * @brief   Test ModelMetadata class for Closed Boundaries (CB)
 * @details
 *
 * Test ModelMetadata class for Closed Boundaries (CB). Check that it correctly reads the
 * partition metadata from a netCDF file and populates the relevant member variables.
 */

#include <doctest/extensions/doctest_mpi.h>
#include <iostream>

#include "ModelMPI.hpp"
#include "ModelMetadata.hpp"

const std::string testFilesDir = TEST_FILES_DIR;
const std::string partitionFilenameCB = testFilesDir + "/halo_cb_test_partition_metadata_3.nc";

namespace Nextsim {

constexpr ModelMetadata::Edge BOTTOM = ModelMetadata::Edge::BOTTOM;
constexpr ModelMetadata::Edge RIGHT = ModelMetadata::Edge::RIGHT;
constexpr ModelMetadata::Edge TOP = ModelMetadata::Edge::TOP;
constexpr ModelMetadata::Edge LEFT = ModelMetadata::Edge::LEFT;

constexpr ModelMetadata::Corner TOP_LEFT = ModelMetadata::Corner::TOP_LEFT;
constexpr ModelMetadata::Corner TOP_RIGHT = ModelMetadata::Corner::TOP_RIGHT;
constexpr ModelMetadata::Corner BOTTOM_RIGHT = ModelMetadata::Corner::BOTTOM_RIGHT;
constexpr ModelMetadata::Corner BOTTOM_LEFT = ModelMetadata::Corner::BOTTOM_LEFT;

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
MPI_TEST_CASE("Test getPartitionMetadata closed boundary", 3)
{
    auto& modelMPI = ModelMPI::getInstance(test_comm);
    auto& meta = ModelMetadata::getInstance(partitionFilenameCB);
    REQUIRE(modelMPI.getComm() == test_comm);
    // this metadata is specific to the non-periodic boundary conditions
    testNonPeriodicBC(test_rank);

    // This metadata is specific to the periodic boundary conditions.
    // They are all zero because the input metadata file `partitionFilenameCB` does not use periodic
    // boundary conditions.
    REQUIRE(meta.neighbourRanksPeriodic[LEFT].size() == 0);
    REQUIRE(meta.neighbourRanksPeriodic[RIGHT].size() == 0);
    REQUIRE(meta.neighbourRanksPeriodic[BOTTOM].size() == 0);
    REQUIRE(meta.neighbourRanksPeriodic[TOP].size() == 0);
    REQUIRE(meta.cornerRanksPeriodic[LEFT].size() == 0);
    REQUIRE(meta.cornerRanksPeriodic[RIGHT].size() == 0);
    REQUIRE(meta.cornerRanksPeriodic[BOTTOM].size() == 0);
    REQUIRE(meta.cornerRanksPeriodic[TOP].size() == 0);
}

}
