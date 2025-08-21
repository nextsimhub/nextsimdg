/*!
 *
 * @author Tom Meltzer <tdm39@cam.ac.uk>
*/

#include <doctest/extensions/doctest_mpi.h>
#include <iostream>

#include "ModelMPI.hpp"
#include "ModelMetadata.hpp"

const std::string testFilesDir = TEST_FILES_DIR;
const std::string partitionFilenameCB = testFilesDir + "/partition_metadata_3_cb.nc";

namespace Nextsim {

constexpr ModelMetadata::Edge BOTTOM = ModelMetadata::Edge::BOTTOM;
constexpr ModelMetadata::Edge RIGHT = ModelMetadata::Edge::RIGHT;
constexpr ModelMetadata::Edge TOP = ModelMetadata::Edge::TOP;
constexpr ModelMetadata::Edge LEFT = ModelMetadata::Edge::LEFT;

typedef std::vector<int> vec;

// these tests are the same for closed boundary conditions (BC) and peridic BC
static void testNonPeriodicBC(int test_rank)
{
    auto& meta = ModelMetadata::getInstance();
    if (test_rank == 0) {
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
    } else if (test_rank == 1) {
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
    } else if (test_rank == 2) {
        REQUIRE(meta.neighbourRanks[LEFT] == vec { 0, 1 });
        REQUIRE(meta.neighbourExtents[LEFT] == vec { 4, 5 });
        REQUIRE(meta.neighbourHaloSend[LEFT] == vec { 7, 7 });
        REQUIRE(meta.neighbourHaloRecv[LEFT] == vec { 15, 19 });
        REQUIRE(meta.neighbourRanks[RIGHT].size() == 0);
        REQUIRE(meta.neighbourRanks[BOTTOM].size() == 0);
        REQUIRE(meta.neighbourRanks[TOP].size() == 0);
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
}

}
