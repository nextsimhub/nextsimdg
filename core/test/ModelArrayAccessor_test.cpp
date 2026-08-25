/*!
 * @author  Tim Spain <timothy.spain@nersc.no>
 * @author  Robert Jendersie <robert.jendersie@ovgu.de>
 */

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "../src/include/ModelArrayAccessor.hpp"
#include "../src/include/ModelArrayStore.hpp"

namespace Nextsim {

class MiniModelComponent {
public:
    static constexpr TextTag H_ICE0 = { "H_ICE0" };
    static constexpr TextTag SW_IN = { "SW_IN" };
    static constexpr TextTag H_ICE = { "H_ICE" };
};

class AtmIn : public MiniModelComponent {
public:
    AtmIn(ModelArrayStore& store)
        : hiceAccessor(store, RW)
        , swinAccessor(store, RW)
    {
    }
    void configure()
    {
        hiceAccessor.getHostRW().reinitialize();
        swinAccessor.getHostRW().reinitialize();
    }
    void setData(const std::vector<FloatType>& values)
    {
        hiceAccessor.getHostRW() = values[0];
        swinAccessor.getHostRW() = values[1];
    }

private:
    ModelArrayAccessor<H_ICE0, RW> hiceAccessor;
    ModelArrayAccessor<SW_IN, RW> swinAccessor;
};

class IceThermo : public MiniModelComponent {
public:
    IceThermo(ModelArrayStore& store)
        : hiceAccessor(store)
    {
    }

    void update(int tStep)
    {
        ModelArray& hice = hiceAccessor.getHostRW();
        hice[0] *= (1. + tStep) / tStep;
    }

private:
    ModelArrayAccessor<H_ICE, RW> hiceAccessor;
};

class IceCalc : public MiniModelComponent {
public:
    IceCalc(ModelArrayStore& store)
        : hiceAccessor(store, RW)
        , hice0Accessor(store)
        , thermo(store)
    {
    }
    void configure() { hiceAccessor.getHostRW().reinitialize(); }
    void update(int tStep)
    {
        ModelArray& hice = hiceAccessor.getHostRW();
        const ModelArray& hice0 = hice0Accessor.getHostRO();
        hice[0] = hice0[0];
        thermo.update(tStep);
    }
    void getData(FloatType& dataOut)
    {
        const ModelArray& hice = hiceAccessor.getHostRO();
        dataOut = hice[0];
    }

private:
    ModelArrayAccessor<H_ICE, RW> hiceAccessor;
    ModelArrayAccessor<H_ICE0> hice0Accessor;

    IceThermo thermo;
};

TEST_SUITE_BEGIN("[ModelArrayAccessor]");
TEST_CASE("Accessing the data")
{
    // create store on stack so that test cases do not influence each other
    ModelArrayStore store;

    AtmIn atmIn(store);
    FloatType hice0 = 0.56;
    FloatType swin = 311;
    ModelArray::setDimensions(ModelArray::Type::H, { 1, 1 });
    atmIn.configure();
    atmIn.setData({ hice0, swin });

    IceCalc iceCalc(store);
    iceCalc.configure();
    int tStep = 40;
    iceCalc.update(tStep);

    FloatType hicef;
    iceCalc.getData(hicef);
    FloatType target = hice0 * (1. + tStep) / tStep;
    REQUIRE(hicef == doctest::Approx(target).epsilon(1e-8));
}

TEST_CASE("(Not) getting write access to a read only field")
{
    ModelArrayStore store;
    ModelArrayAccessor<MiniModelComponent::H_ICE0, RW> hice0SrcAccessor(store, RO);
    HField& hice0Src = hice0SrcAccessor.getHostRW();
    hice0Src.reinitialize();
    hice0Src[0] = 1.0;
    REQUIRE_THROWS_AS(
        (ModelArrayAccessor<MiniModelComponent::H_ICE0, RW>(store)), std::logic_error);

    // inverted initialization order: register after accessor was created
    ModelArrayAccessor<MiniModelComponent::H_ICE, RW> hiceAccessor(store);
    REQUIRE_THROWS_AS(
        (ModelArrayAccessor<MiniModelComponent::H_ICE, RW>(store, RO)), std::logic_error);
}

static const FloatType targetFlux = 320;
// "inline" here prevents warning from -Wsubobject-linkage
// internal linkage which causes the warning would only be a problem if CouplEr was used elsewhere
inline constexpr TextTag sw_in = { "sw_in" };

class CouplEr {
public:
    CouplEr(ModelArrayStore& bs)
        : swFluxAccessor(bs)
    {
    }
    void update()
    {
        ModelArray& swFlux = swFluxAccessor.getHostRW();
        swFlux[0] = targetFlux;
    }

private:
    ModelArrayAccessor<Nextsim::sw_in, RW> swFluxAccessor;
};

class CouplIn : public MiniModelComponent {
public:
    CouplIn(ModelArrayStore& store)
        : hiceAccessor(store, RW)
        , swinAccessor(store, RW)
        , coupler(store)
    {
    }
    void configure()
    {
        hiceAccessor.getHostRW().reinitialize();
        swinAccessor.getHostRW().reinitialize();
    }
    void setData()
    {
        ModelArray& hice = hiceAccessor.getHostRW();
        ModelArray& swin = swinAccessor.getHostRW();
        hice[0] = 0.5;
        swin[0] = 350;
    }
    void update() { coupler.update(); }

private:
    ModelArrayAccessor<H_ICE, RW> hiceAccessor;
    ModelArrayAccessor<SW_IN, RW> swinAccessor;
    CouplEr coupler;
};

TEST_CASE("Accessing the data two ways")
{
    ModelArrayStore store;

    CouplIn couplIn(store);
    ModelArray::setDimensions(ModelArray::Type::H, { 1, 1 });
    couplIn.configure();
    ModelArrayAccessor<sw_in> swinAccessor(store);
    couplIn.setData();

    const ModelArray& swin = swinAccessor.getHostRO();
    REQUIRE(swin[0] != targetFlux);
    couplIn.update();
    REQUIRE(swin[0] == targetFlux);
}

TEST_CASE("Full component access")
{
    ModelArrayStore store;

    const size_t nx = 5;
    const size_t ny = 7;
    const size_t nDG = 3;
    ModelArray::setDimension(ModelArray::Dimension::X, nx);
    ModelArray::setDimension(ModelArray::Dimension::Y, ny);
    ModelArray::setDimension(ModelArray::Dimension::DG, nDG);
    ModelArray dgSrc(ModelArray::Type::DG);
    dgSrc.reinitialize();
    dgSrc = 5.;
    static constexpr TextTag DG_SRC = { "DG_SRC" };
    ModelArrayAccessor<DG_SRC, RW> dgRefAccessor(store, RO, ModelArray::Type::DG);
    ModelArrayAccessor<DG_SRC> dgRefConstAccessor(store);
    const ModelArray& dgRef = dgRefConstAccessor.getHostRO();
    const ModelArray::DataType& eArray = dgRef.data();
    REQUIRE(eArray.rows() == nx * ny);
    REQUIRE(eArray.cols() == nDG);
}

#ifdef USE_KOKKOS
TEST_CASE("Host device data syncing")
{
    Kokkos::initialize(Kokkos::InitializationSettings {});
    // scope to limit the lifetime of the ModelArrayStore
    {
        ModelArrayStore store;

        ModelArray::setDimensions(ModelArray::Type::H, { 2, 3 });

        constexpr int IDX = 5;

        // init on host
        ModelArrayAccessor<MiniModelComponent::H_ICE, RW> hiceSrcAccessor(store, RW);
        {
            ModelArray& hice = hiceSrcAccessor.getHostRW();
            hice.reinitialize();
            hice = 2.0;
            hice[IDX] = 5.0;
        }
        // do work on device
        {
            ModelArrayAccessor<MiniModelComponent::H_ICE, RW> hiceDstAccessor(store);
            const DeviceViewMA& hiceDevice = hiceDstAccessor.getDeviceRW();
            Kokkos::parallel_for(
                "updateDevice", hiceDevice.extent(0),
                KOKKOS_LAMBDA(const DeviceIndex i) { hiceDevice(i, 0) += i == IDX ? -1.0 : 1.0; });
        }
        // check results on host
        {
            const ModelArray& hice = hiceSrcAccessor.getHostRO();
            for (size_t i = 0; i < hice.size(); ++i) {
                if (i == IDX) {
                    CHECK(hice[i] == 4.0);
                } else {
                    CHECK(hice[i] == 3.0);
                }
            }
        }
    }

    Kokkos::finalize();
}
#endif

TEST_SUITE_END();

};
