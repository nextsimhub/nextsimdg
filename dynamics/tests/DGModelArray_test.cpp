/*!
 * @file DGModelArray_test.cpp
 *
 * @date Oct 6, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include "include/DGModelArray.hpp"

#include "include/ModelArrayRef.hpp"
#include "include/ParametricMesh.hpp"

namespace Nextsim {

TEST_CASE("DGVector from ModelArray", "[DGModelArray]")
{
    static const int DG = 3;
    const size_t nx = 32;
    const size_t ny = 32;
    const size_t mx = 100;
    const size_t my = 100;

    ModelArray::setDimensions(ModelArray::Type::DG, { nx, ny });
    ModelArray::setNComponents(ModelArray::Type::DG, DG);

    ModelArray source(ModelArray::Type::DG);
    source.resize();
    // Fill with data
    for (size_t i = 0; i < nx; ++i) {
        for (size_t j = 0; j < ny; j++) {
            for (size_t c = 0; c < DG; ++c) {
                source.components({i, j})[c] = c + my * (j + mx * (i));
            }
        }
    }

    CHECK(source.components({14, 12})[2] == (2 + my * (12 + mx * (14))));

    // Create the ParametricMesh object
    ParametricMesh smash;
    smash.nx = nx;
    smash.ny = ny;
    smash.nnodes = nx * ny;
    smash.nelements = nx * ny;
    smash.vertices.resize(smash.nelements, Eigen::NoChange);
    for (size_t i = 0; i < nx; ++i) {
        for (size_t j = 0; j < ny; ++j) {
            smash.vertices(i * ny + j, 0) = i;
            smash.vertices(i * ny + j, 1) = j;
        }
    }
    DGVector<DG> dest(smash);
    DGModelArray::ma2dg<DG>(source, dest);

    // Did it work?
    size_t targetPoint = ModelArray::indexFromLocation(ModelArray::Type::DG, {14, 12});
    REQUIRE(dest(targetPoint, 2) == source.components({14, 12})[2]);
}

TEST_CASE("ModelArray from DGVector", "[DGModelArray]")
{
    static const int DG = 6;
    const size_t nx = 31;
    const size_t ny = 33;
    const size_t mx = 100;
    const size_t my = 100;

    // Create the ParametricMesh object
    ParametricMesh smash;
    smash.nx = nx;
    smash.ny = ny;
    smash.nnodes = nx * ny;
    smash.nelements = nx * ny;
    smash.vertices.resize(smash.nelements, Eigen::NoChange);
    for (size_t i = 0; i < nx; ++i) {
        for (size_t j = 0; j < ny; ++j) {
            smash.vertices(i * ny + j, 0) = i;
            smash.vertices(i * ny + j, 1) = j;
        }
    }
    DGVector<DG> source(smash);
    // Let the H arrays do the index calculation
    ModelArray::setDimensions(ModelArray::Type::H, { nx, ny });
    for (size_t i = 0; i < nx; ++i) {
        for (size_t j = 0; j < ny; j++) {
            size_t index = ModelArray::indexFromLocation(ModelArray::Type::H, {i, j});
            for (size_t c = 0; c < DG; ++c) {
                source(index, c) = c + my * (j + mx * (i));
            }
        }
    }
    size_t targetPoint = ModelArray::indexFromLocation(ModelArray::Type::H, {14, 12});
    CHECK(source(targetPoint, 4) == 4 + my * (12 + mx *(14)));

    ModelArray::setDimensions(ModelArray::Type::DG, { nx, ny });
    ModelArray::setNComponents(ModelArray::Type::DG, DG);

    ModelArray dest(ModelArray::Type::DG);
    dest.resize();
    DGModelArray::dg2ma(source, dest);

    // Did it work?
    REQUIRE(dest.components({14, 12})[4] == source(targetPoint, 4));
}

TEST_CASE("Test with DG = 1", "[DGModelArray]") // (It would be a silly case to get wrong!)
{
    static const int DG = 1;
    const size_t nx = 23;
    const size_t ny = 27;
    const size_t mx = 100;

    // Create the ParametricMesh object
    ParametricMesh smash;
    smash.nx = nx;
    smash.ny = ny;
    smash.nnodes = nx * ny;
    smash.nelements = nx * ny;
    smash.vertices.resize(smash.nelements, Eigen::NoChange);
    for (size_t i = 0; i < nx; ++i) {
        for (size_t j = 0; j < ny; ++j) {
            smash.vertices(i * ny + j, 0) = i;
            smash.vertices(i * ny + j, 1) = j;
        }
    }
    DGVector<DG> source(smash);
    // Let the H arrays do the index calculation
    ModelArray::setDimensions(ModelArray::Type::H, { nx, ny });
    for (size_t i = 0; i < nx; ++i) {
        for (size_t j = 0; j < ny; j++) {
            size_t index = ModelArray::indexFromLocation(ModelArray::Type::H, {i, j});
            for (size_t c = 0; c < DG; ++c) {
                source(index, c) = j + mx * i;
            }
        }
    }

    size_t targetPoint = ModelArray::indexFromLocation(ModelArray::Type::H, {14, 12});
    CHECK(source(targetPoint) == 12 + mx * 14);

    ModelArray::setDimensions(ModelArray::Type::DG, { nx, ny });
    ModelArray::setNComponents(ModelArray::Type::DG, DG);

    ModelArray dest(ModelArray::Type::DG);
    dest.resize();
    DGModelArray::dg2ma(source, dest);

    // Did it work?
    REQUIRE(dest(14, 12) == source(targetPoint));
    targetPoint = ModelArray::indexFromLocation(ModelArray::Type::H, {12, 14});
    REQUIRE(dest(12, 14) == source(targetPoint));

    // Change the values
    dest *= 2.;
    dest += 1.;

    targetPoint = ModelArray::indexFromLocation(ModelArray::Type::H, {11, 13});
    // Make sure it is different before the copy
    REQUIRE(source(targetPoint) != dest(11, 13));
    DGModelArray::ma2dg<DG>(dest, source);
    // And identical again after
    REQUIRE(source(targetPoint) == dest(11, 13));

}

TEST_CASE("DGVector and HField", "[DGModelArray]")
{
    static const int DG = 3;
    const size_t nx = 29;
    const size_t ny = 27;
    const size_t mx = 100;
    const size_t my = 100;

    ModelArray::setDimensions(ModelArray::Type::H, { nx, ny });
    ModelArray::setDimensions(ModelArray::Type::DG, { nx, ny });
    ModelArray::setNComponents(ModelArray::Type::DG, DG);

    // Create the ParametricMesh object
    ParametricMesh smash;
    smash.nx = nx;
    smash.ny = ny;
    smash.nnodes = nx * ny;
    smash.nelements = nx * ny;
    smash.vertices.resize(smash.nelements, Eigen::NoChange);
    for (size_t i = 0; i < nx; ++i) {
        for (size_t j = 0; j < ny; ++j) {
            smash.vertices(i * ny + j, 0) = i;
            smash.vertices(i * ny + j, 1) = j;
        }
    }

    // THe data comes from the restart file as a Type::DG ModelArray
    ModelArray restart(ModelArray::Type::DG);
    restart.resize();
    // Fill with data
    for (size_t i = 0; i < nx; ++i) {
        for (size_t j = 0; j < ny; j++) {
            for (size_t c = 0; c < DG; ++c) {
                restart.components({i, j})[c] = c + my * (j + mx * (i));
            }
        }
    }

    // This is written to a DGVector stored in the dynamics
    DGVector<DG> dyn(smash);
    DGModelArray::ma2dg<DG>(restart, dyn);

    size_t targetPoint = ModelArray::indexFromLocation(ModelArray::Type::DG, {14, 12});
    REQUIRE(dyn(targetPoint, 2) == restart.components({14, 12})[2]);

    // In a given timestep, the dynamics happens first, changing the DGVector in all components
    dyn(targetPoint, 0) += 0.2983;
    dyn(targetPoint, 1) -= 0.34596;
    dyn(targetPoint, 2) += 0.53244;

    double preservedComponent = dyn(targetPoint, 1);

    // This is the copied out to the thermodynamics array, an HField ModelArray
    ModelArray thermo(ModelArray::Type::H);
    thermo.resize();
    DGModelArray::dg2ma<DG>(dyn, thermo);
    REQUIRE(thermo(14, 12) == dyn(targetPoint, 0));

    // The thermodynamics happens, changing the thermo array
    thermo(14, 12) -= 0.62748;

    // This is then copied back into the dynamics array, only changing the 0 component
    DGModelArray::ma2dg<DG>(thermo, dyn);
    REQUIRE(dyn(targetPoint, 0) == thermo(14, 12));
    REQUIRE(dyn(targetPoint, 1) == preservedComponent);

}

}
