/*!
 * @file ModelArrayRefDebug_test.cpp
 *
 * @date 9 Feb 2023
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "include/ModelArrayRef.hpp"

namespace Nextsim {

class MiniModelComponent {
public:
    enum class ProtectedArray { H_ICE, SW_IN, COUNT };
    enum class SharedArray { H_ICE, COUNT };
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
MARConstBackingStore MiniModelComponent::protectedArrays(
    static_cast<size_t>(ProtectedArray::COUNT));

class IceThermo : public MiniModelComponent {
public:
    IceThermo()
        : hice(MiniModelComponent::getSharedArrays())
    {
    }

    void update(int tStep) { hice[0] *= (1. + tStep) / tStep; }

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
    void configure() { hice.resize(); }
    void update(int tStep)
    {
        hice[0] = hice0[0];
        thermo.update(tStep);
    }
    void getData(double& dataOut) { dataOut = hice[0]; }

private:
    HField hice;
    ModelArrayRef<ProtectedArray::H_ICE, MARConstBackingStore> hice0;

    IceThermo thermo;
};

TEST_SUITE_BEGIN("ModelArrayRefDebug");
TEST_CASE("No registered array")
{
    ModelArray::setDimensions(ModelArray::Type::H, { 1, 1 });

    IceThermo iceThermo;
    REQUIRE_THROWS_AS(iceThermo.update(1), std::invalid_argument);
}

TEST_CASE("Correct access")
{
    ModelArray::setDimensions(ModelArray::Type::H, { 1, 1 });

    IceThermo iceThermo;
    IceCalc iceCalc;
    iceCalc.configure();
    REQUIRE_NOTHROW(iceThermo.update(1));
}
TEST_SUITE_END();
}
