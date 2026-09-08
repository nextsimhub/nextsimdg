/*!
 * @author  Tom Meltzer <tdm39@cam.ac.uk>
 * @brief   Test tripolar partition metadata reading in ModelMetadata
 * @details
 *
 * Verify that ModelMetadata correctly detects tripolar topology from the global "tripolar"
 * attribute written by the domain decomposition tool. This test reads a partition metadata file
 * generated with --tripolar, so the attribute is present. It must run as its own binary because
 * ModelMetadata is a singleton bound to the first partition file it reads; see
 * testModelMetadata_MPI3 (periodic) for the complement of this test.
 */

#include <doctest/extensions/doctest_mpi.h>
#include <iostream>

#include "ModelMPI.hpp"
#include "ModelMetadata.hpp"

const std::string testFilesDir = TEST_FILES_DIR;
const std::string partitionFilename = testFilesDir + "/halo_tp_test_partition_metadata_3.nc";

namespace Nextsim {

TEST_SUITE_BEGIN("ModelMetadata");
MPI_TEST_CASE("Test getPartitionMetadata detects tripolar topology", 3)
{
    auto& modelMPI = ModelMPI::getInstance(test_comm);
    auto& meta = ModelMetadata::getInstance(partitionFilename);

    CHECK(modelMPI.getComm() == test_comm);
    CHECK(meta.usingTripolarTopology());

    switch (test_rank) {
    case 0:
        CHECK(!meta.needsTripolarFold());
        break;
    case 1:
    case 2:
        CHECK(meta.needsTripolarFold());
        break;
    }
}

}
