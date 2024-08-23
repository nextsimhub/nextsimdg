/*!
 * @file ModelMetadata_test.cpp
 *
 * @date 18 Nov 2024
 * @author Tom Meltzer <tdm39@cam.ac.uk>
 */

#include <doctest/extensions/doctest_mpi.h>
#include <iostream>

#include "ModelMetadata.hpp"

const std::string testFilesDir = TEST_FILES_DIR;
const std::string filename = testFilesDir + "/paraGrid_test.nc";
const std::string diagFile = "paraGrid_diag.nc";
const std::string dateString = "2000-01-01T00:00:00Z";
const std::string partitionFilenameCB = testFilesDir + "/partition_metadata_3_cb.nc";
const std::string partitionFilenamePB = testFilesDir + "/partition_metadata_3_pb.nc";

namespace Nextsim {

TEST_SUITE_BEGIN("ModelMetadata");
MPI_TEST_CASE("Test getPartitionMetadata closed boundary", 3)
{
    ModelMetadata metadata(partitionFilenameCB, test_comm);
    REQUIRE(metadata.mpiComm == test_comm);
    // this metadata is specific to the non-periodic boundary conditions
    if (test_rank == 0) {
        REQUIRE(metadata.neighbourRanks[ModelMetadata::LEFT].size() == 0);
        REQUIRE(metadata.neighbourRanks[ModelMetadata::RIGHT] == std::vector<int> { 2 });
        REQUIRE(metadata.neighbourExtents[ModelMetadata::RIGHT] == std::vector<int> { 4 });
        REQUIRE(metadata.neighbourHaloStarts[ModelMetadata::RIGHT] == std::vector<int> { 0 });
        REQUIRE(metadata.neighbourRanks[ModelMetadata::BOTTOM].size() == 0);
        REQUIRE(metadata.neighbourRanks[ModelMetadata::TOP] == std::vector<int> { 1 });
        REQUIRE(metadata.neighbourExtents[ModelMetadata::TOP] == std::vector<int> { 7 });
        REQUIRE(metadata.neighbourHaloStarts[ModelMetadata::TOP] == std::vector<int> { 0 });
    } else if (test_rank == 1) {
        REQUIRE(metadata.neighbourRanks[ModelMetadata::LEFT].size() == 0);
        REQUIRE(metadata.neighbourRanks[ModelMetadata::RIGHT] == std::vector<int> { 2 });
        REQUIRE(metadata.neighbourExtents[ModelMetadata::RIGHT] == std::vector<int> { 5 });
        REQUIRE(metadata.neighbourHaloStarts[ModelMetadata::RIGHT] == std::vector<int> { 12 });
        REQUIRE(metadata.neighbourRanks[ModelMetadata::BOTTOM] == std::vector<int> { 0 });
        REQUIRE(metadata.neighbourExtents[ModelMetadata::BOTTOM] == std::vector<int> { 7 });
        REQUIRE(metadata.neighbourHaloStarts[ModelMetadata::BOTTOM] == std::vector<int> { 21 });
        REQUIRE(metadata.neighbourRanks[ModelMetadata::TOP].size() == 0);
    } else if (test_rank == 2) {
        REQUIRE(metadata.neighbourRanks[ModelMetadata::LEFT] == std::vector<int> { 0, 1 });
        REQUIRE(metadata.neighbourExtents[ModelMetadata::LEFT] == std::vector<int> { 4, 5 });
        REQUIRE(metadata.neighbourHaloStarts[ModelMetadata::LEFT] == std::vector<int> { 6, 6 });
        REQUIRE(metadata.neighbourRanks[ModelMetadata::RIGHT].size() == 0);
        REQUIRE(metadata.neighbourRanks[ModelMetadata::BOTTOM].size() == 0);
        REQUIRE(metadata.neighbourRanks[ModelMetadata::TOP].size() == 0);
    } else {
        std::cerr << "only valid for 3 ranks" << std::endl;
        exit(1);
    }

    // This metadata is specific to the periodic boundary conditions.
    // They are all zero because the input metadata file `partitionFilenameCB` does not use periodic
    // boundary conditions.
    REQUIRE(metadata.neighbourRanksPeriodic[ModelMetadata::LEFT].size() == 0);
    REQUIRE(metadata.neighbourRanksPeriodic[ModelMetadata::RIGHT].size() == 0);
    REQUIRE(metadata.neighbourRanksPeriodic[ModelMetadata::BOTTOM].size() == 0);
    REQUIRE(metadata.neighbourRanksPeriodic[ModelMetadata::TOP].size() == 0);
}

MPI_TEST_CASE("Test getPartitionMetadata periodic boundary", 3)
{
    ModelMetadata metadata(partitionFilenamePB, test_comm);
    REQUIRE(metadata.mpiComm == test_comm);
    // this metadata should be identical to the Closed Boundary version so we check it again
    if (test_rank == 0) {
        REQUIRE(metadata.neighbourRanks[ModelMetadata::LEFT].size() == 0);
        REQUIRE(metadata.neighbourRanks[ModelMetadata::RIGHT] == std::vector<int> { 2 });
        REQUIRE(metadata.neighbourExtents[ModelMetadata::RIGHT] == std::vector<int> { 4 });
        REQUIRE(metadata.neighbourHaloStarts[ModelMetadata::RIGHT] == std::vector<int> { 0 });
        REQUIRE(metadata.neighbourRanks[ModelMetadata::BOTTOM].size() == 0);
        REQUIRE(metadata.neighbourRanks[ModelMetadata::TOP] == std::vector<int> { 1 });
        REQUIRE(metadata.neighbourExtents[ModelMetadata::TOP] == std::vector<int> { 7 });
        REQUIRE(metadata.neighbourHaloStarts[ModelMetadata::TOP] == std::vector<int> { 0 });
    } else if (test_rank == 1) {
        REQUIRE(metadata.neighbourRanks[ModelMetadata::LEFT].size() == 0);
        REQUIRE(metadata.neighbourRanks[ModelMetadata::RIGHT] == std::vector<int> { 2 });
        REQUIRE(metadata.neighbourExtents[ModelMetadata::RIGHT] == std::vector<int> { 5 });
        REQUIRE(metadata.neighbourHaloStarts[ModelMetadata::RIGHT] == std::vector<int> { 12 });
        REQUIRE(metadata.neighbourRanks[ModelMetadata::BOTTOM] == std::vector<int> { 0 });
        REQUIRE(metadata.neighbourExtents[ModelMetadata::BOTTOM] == std::vector<int> { 7 });
        REQUIRE(metadata.neighbourHaloStarts[ModelMetadata::BOTTOM] == std::vector<int> { 21 });
        REQUIRE(metadata.neighbourRanks[ModelMetadata::TOP].size() == 0);
    } else if (test_rank == 2) {
        REQUIRE(metadata.neighbourRanks[ModelMetadata::LEFT] == std::vector<int> { 0, 1 });
        REQUIRE(metadata.neighbourExtents[ModelMetadata::LEFT] == std::vector<int> { 4, 5 });
        REQUIRE(metadata.neighbourHaloStarts[ModelMetadata::LEFT] == std::vector<int> { 6, 6 });
        REQUIRE(metadata.neighbourRanks[ModelMetadata::RIGHT].size() == 0);
        REQUIRE(metadata.neighbourRanks[ModelMetadata::BOTTOM].size() == 0);
        REQUIRE(metadata.neighbourRanks[ModelMetadata::TOP].size() == 0);
    } else {
        std::cerr << "only valid for 3 ranks" << std::endl;
        exit(1);
    }

    // this metadata is specific to the periodic boundary conditions
    if (test_rank == 0) {
        // clang-format off
        REQUIRE(metadata.neighbourRanksPeriodic[ModelMetadata::LEFT] == std::vector<int> { 2 });
        REQUIRE(metadata.neighbourExtentsPeriodic[ModelMetadata::LEFT] == std::vector<int> { 4 });
        REQUIRE(metadata.neighbourHaloStartsPeriodic[ModelMetadata::LEFT] == std::vector<int> { 2 });

        REQUIRE(metadata.neighbourRanksPeriodic[ModelMetadata::RIGHT].size() == 0);

        REQUIRE(metadata.neighbourRanksPeriodic[ModelMetadata::BOTTOM] == std::vector<int> { 1 });
        REQUIRE(metadata.neighbourExtentsPeriodic[ModelMetadata::BOTTOM] == std::vector<int> { 7 });
        REQUIRE(metadata.neighbourHaloStartsPeriodic[ModelMetadata::BOTTOM] == std::vector<int> { 28 });

        REQUIRE(metadata.neighbourRanksPeriodic[ModelMetadata::TOP].size() == 0);
    } else if (test_rank == 1) {
        REQUIRE(metadata.neighbourRanksPeriodic[ModelMetadata::LEFT] == std::vector<int> { 2 });
        REQUIRE(metadata.neighbourExtentsPeriodic[ModelMetadata::LEFT] == std::vector<int> { 5 });
        REQUIRE(metadata.neighbourHaloStartsPeriodic[ModelMetadata::LEFT] == std::vector<int> { 14 });

        REQUIRE(metadata.neighbourRanksPeriodic[ModelMetadata::RIGHT].size() == 0);

        REQUIRE(metadata.neighbourRanksPeriodic[ModelMetadata::BOTTOM].size() == 0);

        REQUIRE(metadata.neighbourRanksPeriodic[ModelMetadata::TOP] == std::vector<int> { 0 });
        REQUIRE(metadata.neighbourExtentsPeriodic[ModelMetadata::TOP] == std::vector<int> { 7 });
        REQUIRE(metadata.neighbourHaloStartsPeriodic[ModelMetadata::TOP] == std::vector<int> { 0 });
    } else if (test_rank == 2) {
        REQUIRE(metadata.neighbourRanksPeriodic[ModelMetadata::LEFT].size() == 0);

        REQUIRE(metadata.neighbourRanksPeriodic[ModelMetadata::RIGHT] == std::vector<int> { 0, 1 });
        REQUIRE(metadata.neighbourExtentsPeriodic[ModelMetadata::RIGHT] == std::vector<int> { 4, 5 });
        REQUIRE(metadata.neighbourHaloStartsPeriodic[ModelMetadata::RIGHT] == std::vector<int> { 0, 0 });

        REQUIRE(metadata.neighbourRanksPeriodic[ModelMetadata::BOTTOM] == std::vector<int> { 2 });
        REQUIRE(metadata.neighbourExtentsPeriodic[ModelMetadata::BOTTOM] == std::vector<int> { 3 });
        REQUIRE(metadata.neighbourHaloStartsPeriodic[ModelMetadata::BOTTOM] == std::vector<int> { 24 });

        REQUIRE(metadata.neighbourRanksPeriodic[ModelMetadata::TOP] == std::vector<int> { 2 });
        REQUIRE(metadata.neighbourExtentsPeriodic[ModelMetadata::TOP] == std::vector<int> { 3 });
        REQUIRE(metadata.neighbourHaloStartsPeriodic[ModelMetadata::TOP] == std::vector<int> { 0 });
        // clang-format on
    } else {
        std::cerr << "only valid for 3 ranks" << std::endl;
        exit(1);
    }
}

}
