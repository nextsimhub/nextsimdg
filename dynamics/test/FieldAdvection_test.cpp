/*!
 * @file FieldAdvection_test.cpp
 *
 * @date Feb 26, 2025
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "dgVectorHolder.hpp"
#include "dgLimit.hpp"
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

const size_t snowL = 6;
const size_t snowB = 2;
const size_t snowR = 14;
const size_t snowT = 18;

class TestData {
public:
    using Slice = ArraySlicer::Slice;
    TestData()
    {
    }

    static void init() {
        ModelArray::setDimension(ModelArray::Dimension::X, nx);
        ModelArray::setDimension(ModelArray::Dimension::Y, ny);
        ModelArray::setDimension(ModelArray::Dimension::XVERTEX, nx + 1);
        ModelArray::setDimension(ModelArray::Dimension::YVERTEX, ny + 1);
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
        static ArraySlicer::Slice iceSlice {{{iceL, iceR}, {iceB, iceT}}};
        return iceSlice;
    }
    static ModelArray getCice()
    {
        ModelArray cice(ModelArray::Type::DG);
        cice = 0.;
        cice.component(0) = iceWithBorders(1.0, 0.5);
        return cice;
    }

    static ModelArray getHice()
    {
        ModelArray hice(ModelArray::Type::DG);
        hice = 0.;
        hice.component(0) = iceWithBorders(1.0, 0.25);
        return hice;
    }
    static ModelArray iceWithBorders(double centralValue, double edgeMul)
    {
        double edgeValue = centralValue * edgeMul;
        double vertexValue = edgeValue * edgeMul;
        ModelArray ice(ModelArray::Type::H);
        ice = 0.;
        ice[getIceAreaSlice()] = centralValue;
        Slice lSlice {{{iceL}, {iceB, iceT}}};
        Slice rSlice {{{iceR}, {iceB, iceT}}};
        Slice tSlice {{{iceL, iceR}, {iceT}}};
        Slice bSlice {{{iceL, iceR}, {iceB}}};
        ice[lSlice] = edgeValue;
        ice[rSlice] = edgeValue;
        ice[tSlice] = edgeValue;
        ice[bSlice] = edgeValue;
        ice(iceL, iceT) = vertexValue;
        ice(iceR, iceT) = vertexValue;
        ice(iceL, iceB) = vertexValue;
        ice(iceR, iceB) = vertexValue;

        return ice;
    }

    static ModelArray getHsnow()
    {
        ModelArray hsnow(ModelArray::Type::DG);
        hsnow = 0.;
//        SLice midSlice =
        hsnow[{{{snowL, snowR}, {snowB, snowT}}}] = 1.0;
        hsnow[{{{snowL, snowR}, {snowB}}}] = 0.5;
        hsnow[{{{snowL, snowR}, {snowT-1}}}] = 0.5;

        return hsnow;
    }

    static ModelArray getCoords()
    {
        double ddeg = 0.1; //˚
        double lat0 = -1;

        ModelArray coords(ModelArray::Type::VERTEX);
        // The coords array uses degrees
        coords.resize();
        for (size_t ii = 0; ii < (nx+1)*(ny+1); ++ii) {
            size_t i = ii % (nx + 1);
            size_t j = ii / (nx + 1);
            coords.components(ii)[0] = ddeg * (i - 0.5);
            coords.components(ii)[1] = ddeg * (j - 0.5) + lat0;
        }
        return coords;
    }
};

TEST_SUITE_BEGIN("Field advection");
TEST_CASE("advection")
{
    TestData::init();
    // Land mask
    ModelArray mask = TestData::getMask();
    REQUIRE(mask(0, 0) == 0);
    REQUIRE(mask(1, 1) == 1);
    // Coordinate array
    ModelArray coords = TestData::getCoords();
    REQUIRE(coords.components({0, 0})[0] == -0.05);
    REQUIRE(coords.components({0, 0})[1] == -1.05);
    REQUIRE(coords.components({1, 0})[0] == 0.05);
    REQUIRE(coords.components({1, 0})[1] == -1.05);
    REQUIRE(coords.components({0, 1})[0] == -0.05);
    REQUIRE(coords.components({0, 1})[1] == -0.95);

    // Ice extents
    ModelArray cice = TestData::getCice();
    REQUIRE(cice(0, 0) == 0);
    REQUIRE(cice(1, 1) == 0);
    REQUIRE(cice(2, 2) == 0.25);
    REQUIRE(cice(3, 3) == 1.0);
    REQUIRE(cice(17, 17) == 1.0);
    REQUIRE(cice(18, 17) == 0.5);
    REQUIRE(cice(19, 19) == 0);
    REQUIRE(cice.components({3, 3})[1] == 0.);
    ModelArray cice0 = cice;
    ModelArray hice = TestData::getHice();
    REQUIRE(hice(0, 0) == 0);
    REQUIRE(hice(1, 1) == 0);
    REQUIRE(hice(2, 2) == 0.0625);
    REQUIRE(hice(3, 3) == 1.0);
    REQUIRE(hice(17, 17) == 1.0);
    REQUIRE(hice(18, 17) == 0.25);
    REQUIRE(hice(19, 19) == 0);
    REQUIRE(cice.components({3, 3})[1] == 0.);
    ModelArray hice0 = hice;

    ModelArray hsnow = TestData::getHsnow();
    REQUIRE(hsnow(0, 0) == 0.);
    REQUIRE(hsnow(snowL-1, snowB) == 0.);
    REQUIRE(hsnow(snowL-1, 10) == 0.);
    REQUIRE(hsnow(snowL, 10) == 1.0);
    REQUIRE(hsnow(snowR - 1, 10) == 1.0);
    REQUIRE(hsnow(snowR, 10) == 0.0);
    ModelArray hsnow0 = hsnow;

    double u0 = 10.;
    // Ice velocities
    ModelArray uIce = TestData::iceWithBorders(u0, 1.0);
    ModelArray vIce = TestData::iceWithBorders(0.0, 1.0);
    REQUIRE(uIce(0, 0) == 0);
    REQUIRE(uIce(1, 1) == 0);
    REQUIRE(uIce(2, 2) == u0);
    REQUIRE(uIce(3, 3) == u0);
    REQUIRE(uIce(17, 17) == u0);
    REQUIRE(uIce(18, 17) == u0);
    REQUIRE(uIce(19, 19) == 0);
    // Assume vIce is fine, since it is zero
    // Forcing fields
    ModelArray uWind(ModelArray::Type::H);
    uWind = u0;
    ModelArray zeroArray(ModelArray::Type::H);
    zeroArray = 0.0;

    ModelArray uOcean(ModelArray::Type::H);
    uOcean = u0;

    // Initialise the kernel mesh
    FreeDriftDynamicsKernel<dg> kernel((DynamicsParameters()));
    kernel.initialise(coords, true, mask);
    // Set the data
    kernel.setData(sshName, zeroArray);
    kernel.setData(uName, uIce);
    kernel.setData(vName, zeroArray);
    kernel.setData(uWindName, uWind);
    kernel.setData(vWindName, zeroArray);
    kernel.setData(uOceanName, uOcean);
    kernel.setData(vOceanName, zeroArray);

    kernel.setDGArray(hiceName, hice);
    kernel.setDGArray(ciceName, cice);

    double deltaT = 1200; // seconds
    TimestepTime tst;
    tst.step = deltaT;
    ModelArray deltaH;
    ModelArray deltaC;
    ModelArray deltaSnow;
    DGVectorHolder<dg> snowHolder(hsnow);
    size_t nt = 50;
    for (size_t i = 0; i < nt; ++i) {
        kernel.prepareAdvection();
        kernel.advectField(tst.step.seconds(), snowHolder, 0.0);
        kernel.update(tst);

        deltaH = hice - hice0;
        deltaC = cice - cice0;
        deltaSnow = hsnow - hsnow0;
    }
    // Check that the DG components are actually changing
    REQUIRE(hice.components({3, 3})[1] != 0.);
    REQUIRE(deltaC(iceL, 10) < 0.);
    REQUIRE(deltaC(iceR, 10) > 0.);
    REQUIRE(deltaSnow(snowL, 10) < 0.);
    REQUIRE(deltaSnow(snowR, 10) > 0.);
}
TEST_SUITE_END();
}
