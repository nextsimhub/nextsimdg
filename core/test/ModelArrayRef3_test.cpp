/*!
 * @file NewModelArrayRef_test.cpp
 *
 * @date Sep 19, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include "../src/include/ModelArrayRef3.hpp"
#include "../src/include/MARBackingStore.hpp"

namespace Nextsim {

class MiniModelComponent {
public:
    static const std::string H_ICE0;
    static const std::string SW_IN;
    static const std::string H_ICE;

    static MARBackingStore& getSharedArrays() { return sharedArrays; }
protected:
    static MARBackingStore sharedArrays;
};

MARBackingStore MiniModelComponent::sharedArrays;
const std::string MiniModelComponent::H_ICE0 = "H_ICE0";
const std::string MiniModelComponent::SW_IN = "SW_IN";
const std::string MiniModelComponent::H_ICE = "H_ICE";

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
    : hice(H_ICE, getSharedArrays())
    {
    }

    void update(int tStep)
    {
        hice[0] *= (1. + tStep) / tStep;
    }
private:
    ModelArrayRef<MARBackingStore, RW> hice;
};

class IceCalc : public MiniModelComponent {
public:
    IceCalc()
    : hice0(H_ICE0, getSharedArrays())
    {
        sharedArrays.registerArray(H_ICE, &hice, RW);
    }
    void configure()
    {
        hice.resize();
    }
    void update(int tStep)
    {
        hice[0] = hice0[0];
        thermo.update(tStep);
    }
    void getData(double& dataOut)
    {
        dataOut = hice[0];
    }

private:
    HField hice;
    ModelArrayRef<MARBackingStore> hice0;

    IceThermo thermo;
};


TEST_CASE("Accessing the data", "[ModelArrayRef]")
{
    AtmIn atmIn;
    double hice0 = 0.56;
    double swin = 311;
    ModelArray::setDimensions(ModelArray::Type::H, {1,1});
    atmIn.configure();
    atmIn.setData({hice0, swin});

    IceCalc iceCalc;
    iceCalc.configure();
    int tStep = 40;
    iceCalc.update(tStep);

    double hicef;
    iceCalc.getData(hicef);
    double target = hice0 * (1. + tStep) / tStep;
    REQUIRE(hicef == Approx(target).epsilon(1e-8));
}

/*
 * Uncommenting this test case should result in a compile time error, as an RO ref should not be writable.
*/
//TEST_CASE("(Not) writing to protected arrays", "[ModelArrayRef]")
//{
//    ModelArrayRef<IMARBackingStore> hice0(MiniModelComponent::H_ICE0, MiniModelComponent::getSharedArrays());
//    hice0[0] = 3.141592;
//}

/*
 * Uncommenting this test should result in a run time error, specifically a segmentation violation,
 * as the H_ICE0 array is never available as a RW array.
 */
//TEST_CASE("(Not) writing to protected arrays", "[ModelArrayRef]")
//{
//    HField hice0Src;
//    hice0Src.resize();
//    hice0Src[0] = 1.0;
//    MiniModelComponent::getSharedArrays().registerArray(MiniModelComponent::H_ICE0, &hice0Src);
//    ModelArrayRef<IMARBackingStore, RW> hice0(MiniModelComponent::H_ICE0, MiniModelComponent::getSharedArrays());
//    REQUIRE(hice0[0] != 3.141592);
//}


static const double targetFlux = 320;

class CouplEr
{
public:
    CouplEr(MARBackingStore& bs)
    : swFlux("sw_in", bs)
    {
    }
    void update() { swFlux[0] = targetFlux; }
private:
    ModelArrayRef<MARBackingStore, RW> swFlux;
};

class CouplIn : public MiniModelComponent
{
public:
    CouplIn()
    : coupler(coupledFields)
    {
        sharedArrays.registerArray(H_ICE, &hice);
        sharedArrays.registerArray(SW_IN, &swin);
        // Set the address of the swin array in the local reference backing store
        coupledFields.registerArray("sw_in", &swin, RW);
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
    void update()
    {
        coupler.update();
    }
    MARBackingStore& bs() { return coupledFields; }
private:
    HField hice;
    HField swin;
    MARBackingStore coupledFields;
    CouplEr coupler;
    };

TEST_CASE("Accessing the data two ways", "[ModelArrayRef]")
{
    CouplIn couplIn;
    ModelArray::setDimensions(ModelArray::Type::H, {1,1});
    couplIn.configure();
    ModelArrayRef<MARBackingStore> swin("sw_in", couplIn.bs());
    couplIn.setData();

    REQUIRE(swin[0] != targetFlux);
    couplIn.update();
    REQUIRE(swin[0] == targetFlux);
}

};
