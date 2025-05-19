/*!
 * @file ModelMetadata_test.cpp
 *
 * @date 19 May 2025
 * @author Tom Meltzer <tdm39@cam.ac.uk>
 */

#include <doctest/extensions/doctest_mpi.h>
#include <iostream>

#include "ModelMetadata.hpp"

const std::string testFilesDir = TEST_FILES_DIR;
const std::string partitionFilenameCB = testFilesDir + "/partition_metadata_3_cb.nc";
const std::string partitionFilenamePB = testFilesDir + "/partition_metadata_3_pb.nc";

namespace Nextsim {

constexpr ModelMetadata::Edge BOTTOM = ModelMetadata::Edge::BOTTOM;
constexpr ModelMetadata::Edge RIGHT = ModelMetadata::Edge::RIGHT;
constexpr ModelMetadata::Edge TOP = ModelMetadata::Edge::TOP;
constexpr ModelMetadata::Edge LEFT = ModelMetadata::Edge::LEFT;

typedef std::vector<int> vec;

// these tests are the same for closed boundary conditions (BC) and peridic BC
static void testNonPeriodicBC(ModelMetadata& meta, int test_rank)
{
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
    ModelMetadata& meta = ModelMetadata::getInstance(partitionFilenameCB, test_comm);
    REQUIRE(meta.mpiComm == test_comm);
    // this metadata is specific to the non-periodic boundary conditions
    testNonPeriodicBC(meta, test_rank);

    // This metadata is specific to the periodic boundary conditions.
    // They are all zero because the input metadata file `partitionFilenameCB` does not use periodic
    // boundary conditions.
    REQUIRE(meta.neighbourRanksPeriodic[LEFT].size() == 0);
    REQUIRE(meta.neighbourRanksPeriodic[RIGHT].size() == 0);
    REQUIRE(meta.neighbourRanksPeriodic[BOTTOM].size() == 0);
    REQUIRE(meta.neighbourRanksPeriodic[TOP].size() == 0);
}

MPI_TEST_CASE("Test getPartitionMetadata periodic boundary", 3)
{
    ModelMetadata& meta = ModelMetadata::getInstance(partitionFilenamePB, test_comm);
    REQUIRE(meta.mpiComm == test_comm);
    // this metadata should be identical to the Closed Boundary version so we check it again
    testNonPeriodicBC(meta, test_rank);

    // this metadata is specific to the periodic boundary conditions
    if (test_rank == 0) {
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
    } else if (test_rank == 1) {
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
    } else if (test_rank == 2) {
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
    }
}

}
