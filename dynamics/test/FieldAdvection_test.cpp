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
    using Slice = ArraySlicer::Slice;
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
        mask[Slice({{0}, {}})] = 0.;
        mask[Slice({{-1}, {}})] = 0.;
        mask[Slice({{}, {0}})] = 0.;
        mask[Slice({{}, {-1}})] = 0.;

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
        cice[getIceAreaSlice()] = 1.;
        Slice lSlice {{{iceL}, {}}};
        Slice rSlice {{{iceR}, {}}};
        Slice tSlice {{{}, {iceT}}};
        Slice bSlice {{{}, {iceB}}};
        cice[lSlice] = 0.5 * cice[lSlice];
        cice[rSlice] = 0.5 * cice[rSlice];
        cice[tSlice] = 0.5 * cice[tSlice];
        cice[bSlice] = 0.5 * cice[bSlice];

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
    ModelArray cice = TestData::getCice();
    REQUIRE(cice(0, 0) == 0);
    REQUIRE(cice(1, 1) == 0);
    REQUIRE(cice(2, 2) == 0.25);
    REQUIRE(cice(3, 3) == 1.0);
    REQUIRE(cice(17, 17) == 1.0);
    REQUIRE(cice(18, 17) == 0.5);
    REQUIRE(cice(19, 19) == 0);

}
TEST_SUITE_END();
}
