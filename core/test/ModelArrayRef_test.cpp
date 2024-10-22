/*!
 * @file ModelArrayRef_test.cpp
 *
 * @date 7 Sep 2023
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
TEST_SUITE_END();

};
