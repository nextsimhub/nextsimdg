/*!
 * @file ModelArrayRef_test.cpp
 *
 * @date 24 Sep 2024
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "../src/include/ModelArrayRef.hpp"
#include "../src/include/ModelArrayReferenceStore.hpp"

namespace Nextsim {

class MiniModelComponent {
public:
    static constexpr TextTag H_ICE0 = { "H_ICE0" };
    static constexpr TextTag SW_IN = { "SW_IN" };
    static constexpr TextTag H_ICE = { "H_ICE" };

    static ModelArrayReferenceStore& getSharedArrays() { return sharedArrays; }

protected:
    static ModelArrayReferenceStore sharedArrays;
};

ModelArrayReferenceStore MiniModelComponent::sharedArrays;

class AtmIn : public MiniModelComponent {
public:
    AtmIn()
    {
        sharedArrays.registerArray(H_ICE0, &hice, RO);
        sharedArrays.registerArray(SW_IN, &swin, RO);
    }
    void configure()
    {
        hice.resize();
        swin.resize();
    }
    void setData(const std::vector<double>& values)
    {
        hice = values[0];
        swin = values[1];
    }

private:
    HField hice;
    HField swin;
};

class IceThermo : public MiniModelComponent {
public:
    IceThermo()
        : hice(getSharedArrays())
    {
    }

    void update(int tStep) { hice[0] *= (1. + tStep) / tStep; }

private:
    ModelArrayRef<H_ICE, RW> hice;
};

class IceCalc : public MiniModelComponent {
public:
    IceCalc()
        : hice0(getSharedArrays())
    {
        sharedArrays.registerArray(H_ICE, &hice, RW);
    }
    void configure() { hice.resize(); }
    void update(int tStep)
    {
        hice[0] = hice0[0];
        thermo.update(tStep);
    }
    void getData(double& dataOut) { dataOut = hice[0]; }

private:
    HField hice;
    ModelArrayRef<H_ICE0> hice0;

    IceThermo thermo;
};

TEST_SUITE_BEGIN("[ModelArrayRef]");
TEST_CASE("Accessing the data")
{
    AtmIn atmIn;
    double hice0 = 0.56;
    double swin = 311;
    ModelArray::setDimensions(ModelArray::Type::H, { 1, 1 });
    atmIn.configure();
    atmIn.setData({ hice0, swin });

    IceCalc iceCalc;
    iceCalc.configure();
    int tStep = 40;
    iceCalc.update(tStep);

    double hicef;
    iceCalc.getData(hicef);
    double target = hice0 * (1. + tStep) / tStep;
    REQUIRE(hicef == doctest::Approx(target).epsilon(1e-8));
}

/*
 * Uncommenting this test case should result in a compile time error, as an RO ref should not be
 * writable.
 */
// TEST_CASE("(Not) writing to protected arrays", "[ModelArrayRef]")
//{
//     ModelArrayRef<H_ICE0> hice0(MiniModelComponent::getSharedArrays());
//     hice0[0] = 3.141592;
// }

/*
 * Uncommenting this test should result in a run time error, specifically a segmentation violation,
 * as the H_ICE0 array is never available as a RW array.
 */
// TEST_CASE("(Not) writing to protected arrays", "[ModelArrayRef]")
//{
//     HField hice0Src;
//     hice0Src.resize();
//     hice0Src[0] = 1.0;
//     MiniModelComponent::getSharedArrays().registerArray(MiniModelComponent::H_ICE0, &hice0Src);
//     ModelArrayRef<H_ICE0, RW> hice0(MiniModelComponent::getSharedArrays());
//     REQUIRE(hice0[0] != 3.141592);
// }

static const double targetFlux = 320;
static constexpr TextTag sw_in = { "sw_in" };

class CouplEr {
public:
    CouplEr(ModelArrayReferenceStore& bs)
        : swFlux(bs)
    {
    }
    void update() { swFlux[0] = targetFlux; }

private:
    ModelArrayRef<sw_in, RW> swFlux;
};

class CouplIn : public MiniModelComponent {
public:
    CouplIn()
        : coupler(coupledFields)
    {
        sharedArrays.registerArray(H_ICE, &hice);
        sharedArrays.registerArray(SW_IN, &swin);
        // Set the address of the swin array in the local reference backing store
        coupledFields.registerArray(sw_in, &swin, RW);
    }
    void configure()
    {
        hice.resize();
        swin.resize();
    }
    void setData()
    {
        hice[0] = 0.5;
        swin[0] = 350;
    }
    void update() { coupler.update(); }
    ModelArrayReferenceStore& bs() { return coupledFields; }

private:
    HField hice;
    HField swin;
    ModelArrayReferenceStore coupledFields;
    CouplEr coupler;
};

TEST_CASE("Accessing the data two ways")
{
    CouplIn couplIn;
    ModelArray::setDimensions(ModelArray::Type::H, { 1, 1 });
    couplIn.configure();
    ModelArrayRef<sw_in> swin(couplIn.bs());
    couplIn.setData();

    REQUIRE(swin[0] != targetFlux);
    couplIn.update();
    REQUIRE(swin[0] == targetFlux);
}

TEST_CASE("Test component 0-only operations")
{
    ModelArray::setDimension(ModelArray::Dimension::DG, 2);
    ModelArray dgSrc(ModelArray::Type::DG);
    dgSrc.resize();
    dgSrc = 5.;
    static constexpr TextTag DG_SRC = { "DG_SRC" };
    MiniModelComponent::getSharedArrays().registerArray(DG_SRC, &dgSrc, RO);
    ModelArrayRef<DG_SRC> dgRef(MiniModelComponent::getSharedArrays());
    ModelArray argument(ModelArray::Type::H);
    argument.resize();
    argument = 3.;
    ModelArray sum = dgRef + argument;
    REQUIRE(sum.getType() == ModelArray::Type::H);
    REQUIRE(sum(0, 0) == 5. + 3.);

    ModelArray difference = dgRef - argument;
    REQUIRE(difference.getType() == ModelArray::Type::H);
    REQUIRE(difference(0, 0) == 5. - 3.);

    ModelArray product = dgRef * argument;
    REQUIRE(product.getType() == ModelArray::Type::H);
    REQUIRE(product(0, 0) == 5. * 3.);

    ModelArray ratio = dgRef / argument;
    REQUIRE(ratio.getType() == ModelArray::Type::H);
    REQUIRE(ratio(0, 0) == 5. / 3.);

    double scalar = 3.;
    ModelArray sumScalar = dgRef + scalar;
    REQUIRE(sumScalar.getType() == ModelArray::Type::H);
    REQUIRE(sumScalar(0, 0) == 5. + 3.);

    ModelArray differenceScalar = dgRef - scalar;
    REQUIRE(differenceScalar.getType() == ModelArray::Type::H);
    REQUIRE(differenceScalar(0, 0) == 5. - 3.);

    ModelArray productScalar = dgRef * scalar;
    REQUIRE(productScalar.getType() == ModelArray::Type::H);
    REQUIRE(productScalar(0, 0) == 5. * 3.);

    ModelArray ratioScalar = dgRef / scalar;
    REQUIRE(ratioScalar.getType() == ModelArray::Type::H);
    REQUIRE(ratioScalar(0, 0) == 5. / 3.);

    static constexpr TextTag RW_SRC = { "RW_SRC" };
    MiniModelComponent::getSharedArrays().registerArray(RW_SRC, &dgSrc, RW);
    ModelArrayRef<RW_SRC, RW> rwRef(MiniModelComponent::getSharedArrays());
    argument = 7.;
    ModelArray sumRW = rwRef + argument;
    REQUIRE(sumRW.getType() == ModelArray::Type::H);
    REQUIRE(sumRW(0, 0) == 5. + 7.);

    ModelArray differenceRW = rwRef - argument;
    REQUIRE(differenceRW.getType() == ModelArray::Type::H);
    REQUIRE(differenceRW(0, 0) == 5. - 7.);

    ModelArray productRW = rwRef * argument;
    REQUIRE(productRW.getType() == ModelArray::Type::H);
    REQUIRE(productRW(0, 0) == 5. * 7.);

    ModelArray ratioRW = rwRef / argument;
    REQUIRE(ratioRW.getType() == ModelArray::Type::H);
    REQUIRE(ratioRW(0, 0) == 5. / 7.);

    scalar = 7.;
    ModelArray sumRWScalar = rwRef + scalar;
    REQUIRE(sumRWScalar.getType() == ModelArray::Type::H);
    REQUIRE(sumRWScalar(0, 0) == 5. + 7.);

    ModelArray differenceRWScalar = rwRef - scalar;
    REQUIRE(differenceRWScalar.getType() == ModelArray::Type::H);
    REQUIRE(differenceRWScalar(0, 0) == 5. - 7.);

    ModelArray productRWScalar = rwRef * scalar;
    REQUIRE(productRWScalar.getType() == ModelArray::Type::H);
    REQUIRE(productRWScalar(0, 0) == 5. * 7.);

    ModelArray ratioRWScalar = rwRef / scalar;
    REQUIRE(ratioRWScalar.getType() == ModelArray::Type::H);
    REQUIRE(ratioRWScalar(0, 0) == 5. / 7.);
}

TEST_CASE("Full component access")
{
    const size_t nx = 5;
    const size_t ny = 7;
    const size_t nDG = 3;
    ModelArray::setDimension(ModelArray::Dimension::X, nx);
    ModelArray::setDimension(ModelArray::Dimension::Y, ny);
    ModelArray::setDimension(ModelArray::Dimension::DG, nDG);
    ModelArray dgSrc(ModelArray::Type::DG);
    dgSrc.resize();
    dgSrc = 5.;
    static constexpr TextTag DG_SRC = { "DG_SRC" };
    MiniModelComponent::getSharedArrays().registerArray(DG_SRC, &dgSrc, RO);
    ModelArrayRef<DG_SRC> dgRef(MiniModelComponent::getSharedArrays());
    const ModelArray::DataType& eArray = dgRef.allComponents();
    REQUIRE(eArray.rows() == nx * ny);
    REQUIRE(eArray.cols() == nDG);

}
TEST_SUITE_END();

};
