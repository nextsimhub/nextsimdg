/*!
 * @file NewModelArrayRef_test.cpp
 *
 * @date Sep 19, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include "include/ModelArrayRef2.hpp"

#include <iostream>

namespace Nextsim {

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
    static MARConstBackingStore& getProtectedArrays() { return protectedArrays; }
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
    : hice(MiniModelComponent::sharedArrays)
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
    : hice0(MiniModelComponent::protectedArrays)
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
    ModelArray::setDimensions(ModelArray::Type::H, {1});
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

}


