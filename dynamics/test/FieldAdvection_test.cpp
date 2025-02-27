/*!
 * @file FieldAdvection_test.cpp
 *
 * @date Feb 26, 2025
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "FreeDriftDynamicsKernel.hpp"
#include "include/ModelArray.hpp"
#include "include/ModelArraySlice.hpp"
#include "ParametricMesh.hpp"
#include "include/constants.hpp"

namespace Nextsim {

const size_t nx = 50;
const size_t ny = 21;
const size_t dg = 3;

const size_t iceL = 2;
const size_t iceB = 2;
const size_t iceR = 18;
const size_t iceT = 18;

class TestData {
public:
    TestData()
    {
    }

    static void init() {
        ModelArray::setDimension(ModelArray::Dimension::X, nx);
        ModelArray::setDimension(ModelArray::Dimension::Y, ny);
        ModelArray::setDimension(ModelArray::Dimension::DG, dg);
    }

    static ModelArray getMask()
    {
        ModelArray mask(ModelArray::Type::H);
        // Ocean everywhere
        mask = 1.;
        // Apart from the edges
        ModelArraySlice(mask, {{{0}, {}}}) = 0.;
        ModelArraySlice(mask, {{{-1}, {}}}) = 0.;
        ModelArraySlice(mask, {{{}, {0}}}) = 0.;
        ModelArraySlice(mask, {{{}, {-1}}}) = 0.;

        return mask;
    }

    static ArraySlicer::Slice getIceAreaSlice()
    {
        static ArraySlicer::Slice iceSlice {{{iceL, iceR}, {iceL, iceR}}};
        return iceSlice;
    }
    static ModelArray getCice()
    {
        ModelArray cice(ModelArray::Type::H);
        cice = 0.;
        ModelArraySlice(cice, getIceAreaSlice()) = 1.;
        ModelArraySlice(cice, {{{iceL}, {}}}) *= 0.5;
        ModelArraySlice(cice, {{{iceR - 1}, {}}}) *= 0.5;
        ModelArraySlice(cice, {{{}, {iceB}}}) *= 0.5;
        ModelArraySlice(cice, {{{}, {iceT - 1}}}) *= 0.5;

        return cice;
    }
};

class TestMesh {
    ParametricMesh getMesh()
    {
        ParametricMesh smesh(SPHERICAL);
        smesh.nx = nx;
        smesh.ny = ny;
        smesh.nelements = nx * ny;
        smesh.nnodes = (nx + 1) * (ny + 1);

        smesh.vertices.resize(smesh.nnodes, 2);
        double ddeg = 0.1; //˚
        double d = PhysicalConstants::deg2rad * ddeg;
        double lon0 = -1;
        for (size_t ii = 0; ii < smesh.nnodes; ++ii) {
            size_t i = ii % (nx + 1);
            size_t j = ii / (nx + 1);
            smesh.vertices(ii, 0) = d * (i - 0.5);
            smesh.vertices(ii, 1) = d * (j - 0.5) + lon0;
        }
        return smesh;
    }
};

TEST_SUITE_BEGIN("Field advection");
TEST_CASE("advection")
{
    TestData::init();
    ModelArray mask = TestData::getMask();
    REQUIRE(mask(0, 0) == 0);
    REQUIRE(mask(1, 1) == 1);
}
TEST_SUITE_END();
}
