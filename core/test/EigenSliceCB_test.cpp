/*!
 * @file ModelMetadata_test.cpp
 *
 * @date 16 Jul 2025
 * @author Tom Meltzer <tdm39@cam.ac.uk>
 */

#include <doctest/extensions/doctest_mpi.h>
#include <iostream>

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
MPI_TEST_CASE("test EigenSlice for HField", 3)
{

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

    // create example 2D field on each process
    auto source = ModelArray::HField();
    auto target = ModelArray::HField();
    source.resize();
    target.resize();
    target = -1;

    initializeData(source, test_rank);

    // check we can copy the middle block from one to the other
    Slice sl { { { 1, -1 }, { 1, -1 } } };
    EigenSlice<ModelArray::DataType> sourceSlice(source, sl);
    EigenSlice<ModelArray::DataType> targetSlice(target, sl);

    sourceSlice.copyToSlicedBuffer(targetSlice);
    if (debug) {
        std::cout << "HField" << std::endl;
        targetSlice.print();
    }

    REQUIRE(target(2, 3) == source(2, 3));
    REQUIRE(target(0, 3) == -1);
}

MPI_TEST_CASE("test EigenSlice for DGField", 3)
{
    auto& modelMPI = ModelMPI::getInstance(test_comm);
    auto& metadata = ModelMetadata::getInstance();

    const size_t nx = metadata.getGlobalExtentX();
    const size_t ny = metadata.getGlobalExtentY();
    const size_t localNx = metadata.getLocalExtentX() + 2 * Halo::haloWidth;
    const size_t localNy = metadata.getLocalExtentY() + 2 * Halo::haloWidth;
    const size_t offsetX = metadata.getLocalCornerX();
    const size_t offsetY = metadata.getLocalCornerY();

    ModelArray::setDimension(ModelArray::Dimension::X, nx, localNx, offsetX);
    ModelArray::setDimension(ModelArray::Dimension::Y, ny, localNy, offsetY);
    ModelArray::setDimension(ModelArray::Dimension::DG, DG, DG, 0);

    // create example 2D field on each process
    auto source = ModelArray::DGField();
    auto target = ModelArray::DGField();
    source.resize();
    target.resize();
    target = -1;

    initializeData(source, test_rank);

    // check we can copy the middle block from one to the other
    Slice sl { { { 1, -1 }, { 1, -1 } } };
    EigenSlice<ModelArray::DataType> sourceSlice(source, sl);
    EigenSlice<ModelArray::DataType> targetSlice(target, sl);

    sourceSlice.copyToSlicedBuffer(targetSlice);
    if (debug) {
        std::cout << "DGField" << std::endl;
        targetSlice.print();
    }

    REQUIRE(target.components({ 2, 3 })[1] == source.components({ 2, 3 })[1]);
    REQUIRE(target.components({ 0, 3 })[1] == -1);
}

MPI_TEST_CASE("test EigenSlice for DGVector", 3)
{
    auto& modelMPI = ModelMPI::getInstance(test_comm);
    auto& metadata = ModelMetadata::getInstance();

    const size_t nx = metadata.getGlobalExtentX();
    const size_t ny = metadata.getGlobalExtentY();
    const size_t localNx = metadata.getLocalExtentX() + 2 * Halo::haloWidth;
    const size_t localNy = metadata.getLocalExtentY() + 2 * Halo::haloWidth;
    const size_t offsetX = metadata.getLocalCornerX();
    const size_t offsetY = metadata.getLocalCornerY();

    ModelArray::setDimension(ModelArray::Dimension::X, nx, localNx, offsetX);
    ModelArray::setDimension(ModelArray::Dimension::Y, ny, localNy, offsetY);
    ModelArray::setDimension(ModelArray::Dimension::DG, DG, DG, 0);

    auto fake_data = ModelArray::DGField();
    fake_data.resize();
    initializeData(fake_data, test_rank);

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
    DGVector<DG> source(smesh);
    DGModelArray::ma2dg(fake_data, source);
    DGVector<DG> target(smesh);
    fake_data = -1.;
    DGModelArray::ma2dg(fake_data, target);

    // check we can copy the middle block from one to the other
    Slice sl { { { 1, -1 }, { 1, -1 } } };
    EigenSlice<DGVector<DG>::EigenDGVector> sourceSlice(source, sl, smesh);
    EigenSlice<DGVector<DG>::EigenDGVector> targetSlice(target, sl, smesh);

    sourceSlice.copyToSlicedBuffer(targetSlice);
    if (debug) {
        std::cout << "DGVector" << std::endl;
        targetSlice.print();
    }

    REQUIRE(target(localNx + 1, 1) == source(localNx + 1, 1));
    REQUIRE(target(0, 0) == -1);
}
}
