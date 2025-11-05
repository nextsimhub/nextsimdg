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
        hiceAccessor.getHostRW().resize();
        swinAccessor.getHostRW().resize();
    }
    void setData(const std::vector<double>& values)
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
    void configure() { hiceAccessor.getHostRW().resize(); }
    void update(int tStep)
    {
        ModelArray& hice = hiceAccessor.getHostRW();
        const ModelArray& hice0 = hice0Accessor.getHostRO();
        hice[0] = hice0[0];
        thermo.update(tStep);
    }
    void getData(double& dataOut)
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
    double hice0 = 0.56;
    double swin = 311;
    ModelArray::setDimensions(ModelArray::Type::H, { 1, 1 });
    atmIn.configure();
    atmIn.setData({ hice0, swin });

    IceCalc iceCalc(store);
    iceCalc.configure();
    int tStep = 40;
    iceCalc.update(tStep);

    double hicef;
    iceCalc.getData(hicef);
    double target = hice0 * (1. + tStep) / tStep;
    REQUIRE(hicef == doctest::Approx(target).epsilon(1e-8));
}

TEST_CASE("(Not) getting write access to a read only field")
{
    ModelArrayStore store;
    ModelArrayAccessor<MiniModelComponent::H_ICE0, RW> hice0SrcAccessor(store, RO);
    HField& hice0Src = hice0SrcAccessor.getHostRW();
    hice0Src.resize();
    hice0Src[0] = 1.0;
    REQUIRE_THROWS_AS(
        (ModelArrayAccessor<MiniModelComponent::H_ICE0, RW>(store)), std::logic_error);

    // inverted initialization order: register after accessor was created
    ModelArrayAccessor<MiniModelComponent::H_ICE, RW> hiceAccessor(store);
    REQUIRE_THROWS_AS(
        (ModelArrayAccessor<MiniModelComponent::H_ICE, RW>(store, RO)), std::logic_error);
}

static const double targetFlux = 320;
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
        hiceAccessor.getHostRW().resize();
        swinAccessor.getHostRW().resize();
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

/*
 * This functionality is specific to ModelArrayRef and no longer needed.
 */
// TEST_CASE("Test component 0-only operations")
// {
//     ModelArrayStore store;
//
//     static constexpr TextTag DG_SRC = { "DG_SRC" };
//     ModelArray::setDimension(ModelArray::Dimension::DG, 2);
//     ModelArrayAccessor<DG_SRC, RW> dgSrcAccessor(store, RO, ModelArray::Type::DG);
//     ModelArray& dgSrc = dgSrcAccessor.getHostRW();
//     dgSrc.resize();
//     dgSrc = 5.;
//
//     ModelArrayAccessor<DG_SRC> dgRefAccessor(store);
//     const ModelArray& dgRef = dgRefAccessor.getHostRO();
//     ModelArray argument(ModelArray::Type::H);
//     argument.resize();
//     argument = 3.;
//     ModelArray sum = dgRef + argument;
//     REQUIRE(sum.getType() == ModelArray::Type::H);
//     REQUIRE(sum(0, 0) == 5. + 3.);
//
//     ModelArray difference = dgRef - argument;
//     REQUIRE(difference.getType() == ModelArray::Type::H);
//     REQUIRE(difference(0, 0) == 5. - 3.);
//
//     ModelArray product = dgRef * argument;
//     REQUIRE(product.getType() == ModelArray::Type::H);
//     REQUIRE(product(0, 0) == 5. * 3.);
//
//     ModelArray ratio = dgRef / argument;
//     REQUIRE(ratio.getType() == ModelArray::Type::H);
//     REQUIRE(ratio(0, 0) == 5. / 3.);
//
//     double scalar = 3.;
//     ModelArray sumScalar = dgRef + scalar;
//     REQUIRE(sumScalar.getType() == ModelArray::Type::H);
//     REQUIRE(sumScalar(0, 0) == 5. + 3.);
//
//     ModelArray differenceScalar = dgRef - scalar;
//     REQUIRE(differenceScalar.getType() == ModelArray::Type::H);
//     REQUIRE(differenceScalar(0, 0) == 5. - 3.);
//
//     ModelArray productScalar = dgRef * scalar;
//     REQUIRE(productScalar.getType() == ModelArray::Type::H);
//     REQUIRE(productScalar(0, 0) == 5. * 3.);
//
//     ModelArray ratioScalar = dgRef / scalar;
//     REQUIRE(ratioScalar.getType() == ModelArray::Type::H);
//     REQUIRE(ratioScalar(0, 0) == 5. / 3.);
//
//     static constexpr TextTag RW_SRC = { "RW_SRC" };
//     ModelArrayAccessor<RW_SRC, RW> rwRefAccessor(store, RW);
//     const ModelArray& rwRef = rwRefAccessor.getHostRW();
//     argument = 7.;
//     ModelArray sumRW = rwRef + argument;
//     REQUIRE(sumRW.getType() == ModelArray::Type::H);
//     REQUIRE(sumRW(0, 0) == 5. + 7.);
//
//     ModelArray differenceRW = rwRef - argument;
//     REQUIRE(differenceRW.getType() == ModelArray::Type::H);
//     REQUIRE(differenceRW(0, 0) == 5. - 7.);
//
//     ModelArray productRW = rwRef * argument;
//     REQUIRE(productRW.getType() == ModelArray::Type::H);
//     REQUIRE(productRW(0, 0) == 5. * 7.);
//
//     ModelArray ratioRW = rwRef / argument;
//     REQUIRE(ratioRW.getType() == ModelArray::Type::H);
//     REQUIRE(ratioRW(0, 0) == 5. / 7.);
//
//     scalar = 7.;
//     ModelArray sumRWScalar = rwRef + scalar;
//     REQUIRE(sumRWScalar.getType() == ModelArray::Type::H);
//     REQUIRE(sumRWScalar(0, 0) == 5. + 7.);
//
//     ModelArray differenceRWScalar = rwRef - scalar;
//     REQUIRE(differenceRWScalar.getType() == ModelArray::Type::H);
//     REQUIRE(differenceRWScalar(0, 0) == 5. - 7.);
//
//     ModelArray productRWScalar = rwRef * scalar;
//     REQUIRE(productRWScalar.getType() == ModelArray::Type::H);
//     REQUIRE(productRWScalar(0, 0) == 5. * 7.);
//
//     ModelArray ratioRWScalar = rwRef / scalar;
//     REQUIRE(ratioRWScalar.getType() == ModelArray::Type::H);
//     REQUIRE(ratioRWScalar(0, 0) == 5. / 7.);
// }

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
    dgSrc.resize();
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
            hice.resize();
            hice = 2.0;
            hice[IDX] = 5.0;
        }
        // do work on device
        ModelArrayAccessor<MiniModelComponent::H_ICE, RW> hiceDstAccessor(store);
        {
            const DeviceView& hiceDevice = hiceDstAccessor.getDeviceRW();
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
