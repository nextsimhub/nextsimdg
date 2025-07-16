/*!
 * @file ModelMetadata_test.cpp
 *
 * @date 16 Jul 2025
 * @author Tom Meltzer <tdm39@cam.ac.uk>
 */

#include <doctest/extensions/doctest_mpi.h>

#include "ModelMPI.hpp"
#include "ModelMetadata.hpp"
#include "include/DGModelArray.hpp"
#include "include/EigenSlice.hpp"
#include "include/Halo.hpp"
#include "include/ParametricMesh.hpp"

namespace Nextsim {

const std::string testFilesDir = TEST_FILES_DIR;
const std::string file = testFilesDir + "/partition_metadata_3_cb.nc";

static const int DG = 3;
static const bool debug = false;

using Slice = ArraySlicer::Slice;

void initializeData(ModelArray& ma, int test_rank)
{
    auto& metadata = ModelMetadata::getInstance();

    const size_t localNx = metadata.getLocalExtentX() + 2 * Halo::haloWidth;
    const size_t localNy = metadata.getLocalExtentY() + 2 * Halo::haloWidth;
    const size_t offsetX = metadata.getLocalCornerX();
    const size_t offsetY = metadata.getLocalCornerY();

    auto dof = ma.nComponents();

    // initialize with mock data
    for (size_t j = 0; j < localNy; ++j) {
        for (size_t i = 0; i < localNx; ++i) {
            for (size_t k = 0; k < dof; ++k) {
                ma.components({ i, j })[k]
                    = (test_rank + 1) * 100 + (i + offsetX) * 10 + (j + offsetY) + 1000 * (k + 1);
            }
        }
    }
}

TEST_SUITE_BEGIN("EigenSlice tests");
MPI_TEST_CASE("test EigenSlice on a 3 proc grid", 3)
{
    // test for a 3 proc grid
    // ┌─┬─┐
    // │1│ │
    // ├─┤2│
    // │0│ │
    // └─┴─┘ (proc id)

    auto& modelMPI = ModelMPI::getInstance(test_comm);
    auto& metadata = ModelMetadata::getInstance(file);

    const size_t nx = metadata.getGlobalExtentX();
    const size_t ny = metadata.getGlobalExtentY();
    const size_t localNx = metadata.getLocalExtentX() + 2 * Halo::haloWidth;
    const size_t localNy = metadata.getLocalExtentY() + 2 * Halo::haloWidth;
    const size_t offsetX = metadata.getLocalCornerX();
    const size_t offsetY = metadata.getLocalCornerY();

    ModelArray::setDimension(ModelArray::Dimension::X, nx, localNx, offsetX);
    ModelArray::setDimension(ModelArray::Dimension::Y, ny, localNy, offsetY);
    ModelArray::setDimension(ModelArray::Dimension::DG, DG, DG, 0);

    // HField tests
    {
        // create example 2D field on each process
        auto source = ModelArray::HField();
        auto target = ModelArray::HField();
        source.resize();
        target.resize();
        target = -1;

        initializeData(source, test_rank);

        // check we can copy the middle block from one to the other
        Slice sl { { { 1, -1 }, { 1, -1 }, {} } };
        EigenSlice<ModelArray::DataType> sourceSlice(source, sl);
        EigenSlice<ModelArray::DataType> targetSlice(target, sl);

        // sourceSlice.print();
        // targetSlice.print();
        sourceSlice.copyToSlicedBuffer(targetSlice);
        // targetSlice.print();

        REQUIRE(target(2, 3) == source(2, 3));
        REQUIRE(target(0, 3) == -1);
    }

    // DGField tests
    {
        // create example 2D field on each process
        auto source = ModelArray::DGField();
        auto target = ModelArray::DGField();
        source.resize();
        target.resize();
        target = -1;

        initializeData(source, test_rank);

        // check we can copy the middle block from one to the other
        Slice sl { { { 1, -1 }, { 1, -1 }, {} } };
        EigenSlice<ModelArray::DataType> sourceSlice(source, sl);
        EigenSlice<ModelArray::DataType> targetSlice(target, sl);

        sourceSlice.print();
        targetSlice.print();
        sourceSlice.copyToSlicedBuffer(targetSlice);
        targetSlice.print();

        REQUIRE(target.components({ 2, 3 })[1] == source.components({ 2, 3 })[1]);
        REQUIRE(target.components({ 0, 3 })[1] == -1);
    }

    {
        auto source = ModelArray::DGField();
        source.resize();
        initializeData(source, test_rank);

        // check we can copy the middle block from one to the other
        Slice sl { { { 1, -1 }, { 1, -1 }, { 0, -1 } } };
        EigenSlice<ModelArray::DataType> sourceSlice(source, sl);
        sourceSlice.print();

        ParametricMesh smesh(CARTESIAN);
        smesh.nx = localNx;
        smesh.ny = localNy;
        smesh.nnodes = localNx * localNy;
        smesh.nelements = localNx * localNy;
        smesh.vertices.resize(smesh.nelements, Eigen::NoChange);
        for (size_t i = 0; i < localNx; ++i) {
            for (size_t j = 0; j < localNy; ++j) {
                smesh.vertices(i * localNy + j, 0) = i;
                smesh.vertices(i * localNy + j, 1) = j;
            }
        }
        DGVector<DG> target(smesh);
        target.setZero();

        EigenSlice<DGVector<DG>::EigenDGVector> targetSlice(target, sl, smesh);

        auto ele = target(1, 1);
        ele = target(0, 0);

        targetSlice.print();
        targetSlice.copyFromSlicedBuffer(sourceSlice);
        targetSlice.print();

        ele = target(1, 1);
        ele = target(0, 0);

        REQUIRE(target(2, 3) == source.components({ 2, 3 })[1]);
        REQUIRE(target(0, 0) == -1);
    }

    // for (size_t j = 0; j < localNy; j++) {
    //     for (size_t i = 0; i < localNx; i++) {
    //         REQUIRE(testData(i, j) == mockData(i + j * localNx));
    //     }
    // }
}
}
