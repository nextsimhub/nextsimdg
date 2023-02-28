/*!
 * @file NewModelArrayRef_test.cpp
 *
 * @date Sep 19, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include "../src/include/ModelArrayRef3.hpp"

#include <iostream>

namespace Nextsim {

typedef std::vector<ModelArrayReference> MARBackingStore;
typedef std::vector<ModelArrayConstReference> MARConstBackingStore;

class MiniModelComponent {
public:
    enum class ProtectedArray {
        H_ICE,
        SW_IN,
        COUNT
    };
    enum class SharedArray {
        H_ICE,
        COUNT
    };
    static void registerSharedArray(SharedArray type, ModelArray* p)
    {
        sharedArrays[static_cast<size_t>(type)] = p;
    }
    static void registerProtectedArray(ProtectedArray type, ModelArray* p)
    {
        protectedArrays[static_cast<size_t>(type)] = p;
    }
    static const MARConstBackingStore& getProtectedArrays() { return protectedArrays; }
    static const MARBackingStore& getSharedArrays() { return sharedArrays; }
protected:
    static MARBackingStore sharedArrays;
    static MARConstBackingStore protectedArrays;
};

MARBackingStore MiniModelComponent::sharedArrays(static_cast<size_t>(SharedArray::COUNT));
MARConstBackingStore MiniModelComponent::protectedArrays(static_cast<size_t>(ProtectedArray::COUNT));

class AtmIn : public MiniModelComponent {
public:
    AtmIn()
    {
        registerProtectedArray(ProtectedArray::H_ICE, &hice);
        registerProtectedArray(ProtectedArray::SW_IN, &swin);
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
    : hice(MiniModelComponent::getSharedArrays())
    {
    }

    void update(int tStep)
    {
        hice[0] *= (1. + tStep) / tStep;
    }
private:
    ModelArrayRef<SharedArray::H_ICE, MARBackingStore, RW> hice;
};

class IceCalc : public MiniModelComponent {
public:
    IceCalc()
    : hice0(MiniModelComponent::getProtectedArrays())
    {
        registerSharedArray(SharedArray::H_ICE, &hice);
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
    ModelArrayRef<ProtectedArray::H_ICE, MARConstBackingStore> hice0;

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
 * Uncommenting this test case should result in a compile time error, as a ConstBackingStore backed array ref should not be writable.

TEST_CASE("(Not) writing to protected arrays", "[ModelArrayRef]") {
    ModelArrayRef<MiniModelComponent::ProtectedArray::H_ICE, MARConstBackingStore> hice0(MiniModelComponent::getProtectedArrays());
    hice0[0] = 3.141592;
}
*/
enum class couplFields {
    SWIN,
    COUNT
};

static const double targetFlux = 320;

class CouplEr
{
public:
    CouplEr(MARBackingStore& bs)
    : swFlux(bs)
    {
    }
    void update() { swFlux[0] = targetFlux; }
private:
    ModelArrayRef<couplFields::SWIN, MARBackingStore, RW> swFlux;
};

class CouplIn : public MiniModelComponent
{
public:
    CouplIn()
    : coupledFields(static_cast<size_t>(couplFields::COUNT))
    , coupler(coupledFields)
    {
        registerProtectedArray(ProtectedArray::H_ICE, &hice);
        registerProtectedArray(ProtectedArray::SW_IN, &swin);
        // Set the address of the swin array in the local reference backing store
        coupledFields[static_cast<size_t>(couplFields::SWIN)] = &swin;
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
    const MARBackingStore& bs() { return coupledFields; }
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
    ModelArrayRef<couplFields::SWIN, MARBackingStore> swin(couplIn.bs());
    couplIn.setData();

    REQUIRE(swin[0] != targetFlux);
    couplIn.update();
    REQUIRE(swin[0] == targetFlux);
}
};
