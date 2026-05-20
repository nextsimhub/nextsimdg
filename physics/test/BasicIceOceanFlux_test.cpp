/*!
 * @author  Tim Spain <timothy.spain@nersc.no>
 */

// Does this class need testing? Not really, but it got removed from
// FiniteElementFluxes_test and I thought the tests should continue to exist
// somewhere

#include <doctest/doctest.h>

#include "include/BasicIceOceanHeatFlux.hpp"

#include "include/IFreezingPoint.hpp"
#include "include/ModelArray.hpp"
#include "include/NextsimModule.hpp"
#include "include/Time.hpp"
#include "include/UnescoFreezing.hpp"
#include "include/UniformOcean.hpp"
#include "include/constants.hpp"

namespace Nextsim {

TEST_SUITE_BEGIN("BasicIceOceanHeatFlux");
TEST_CASE("Melting conditions")
{
    ModelArray::setDimensions(ModelArray::Type::H, { 1, 1 });

    Module::setImplementation<IFreezingPoint>("Nextsim::UnescoFreezing");
    UniformOcean ocnBdy;
    ocnBdy.setSST(-1.).setSSS(32.).setMLD(10.25).setU(0).setV(0);
    ocnBdy.setData(ModelState().data);

    class ProgData : public ModelComponent {
    public:
        ProgData()
            : hiceAccessor(getStore(), RW)
            , ciceAccessor(getStore(), RW)
            , hsnowAccessor(getStore(), RW)
        {
        }
        std::string getName() const override { return "ProgData"; }

        void setData(const ModelState::DataMap&) override
        {
            noLandMask();
            ciceAccessor.getHostRW()[0] = 0.5;
            hiceAccessor.getHostRW()[0] = 0.1; // Here we are using the cell-averaged thicknesses
            hsnowAccessor.getHostRW()[0] = 0.01;
        }

        ModelArrayAccessor<Shared::H_ICE_DG, RW> hiceAccessor;
        ModelArrayAccessor<Shared::C_ICE_DG, RW> ciceAccessor;
        ModelArrayAccessor<Shared::H_SNOW_DG, RW> hsnowAccessor;
    } iceState;
    iceState.setData(ModelState().data);

    ModelArrayAccessor<Shared::Q_IO, RW> qioAccessor(ModelComponent::getStore());
    qioAccessor.getHostRW().reinitialize();
    TimestepTime tst = { TimePoint("2000-001"), Duration("P0-0T0:10:0") };
    BasicIceOceanHeatFlux biohf;
    biohf.update(tst);

    constexpr FloatType prec = 1e-5;
    REQUIRE(qioAccessor.getHostRO()[0] == doctest::Approx(53717.8).epsilon(prec));
}

TEST_CASE("Freezing conditions")
{
    ModelArray::setDimensions(ModelArray::Type::H, { 1, 1 });

    Module::setImplementation<IFreezingPoint>("Nextsim::UnescoFreezing");
    UniformOcean ocnBdy;
    ocnBdy.setSST(-1.75).setSSS(32.).setMLD(10.25).setU(0).setV(0);
    ocnBdy.setData(ModelState().data);

    class ProgData : public ModelComponent {
    public:
        ProgData() // RO before the Accessor refactor but these fields have to be RW
            : hiceAccessor(getStore(), RW) // RO
            , ciceAccessor(getStore(), RW) // RO
            , hsnowAccessor(getStore(), RW) // RO
        {
        }
        std::string getName() const override { return "ProgData"; }

        void setData(const ModelState::DataMap&) override
        {
            noLandMask();
            ciceAccessor.getHostRW()[0] = 0.5;
            hiceAccessor.getHostRW()[0] = 0.1; // Here we are using the cell-averaged thicknesses
            hsnowAccessor.getHostRW()[0] = 0.01;
        }

        ModelArrayAccessor<Shared::H_ICE_DG, RW> hiceAccessor;
        ModelArrayAccessor<Shared::C_ICE_DG, RW> ciceAccessor;
        ModelArrayAccessor<Shared::H_SNOW_DG, RW> hsnowAccessor;
    } iceState;
    iceState.setData(ModelState().data);

    ModelArrayAccessor<Shared::Q_IO, RW> qioAccessor(ModelComponent::getStore());
    qioAccessor.getHostRW().reinitialize();

    TimestepTime tst = { TimePoint("2000-001"), Duration("P0-0T0:10:0") };
    BasicIceOceanHeatFlux biohf;
    biohf.update(tst);

    constexpr FloatType prec = 1e-5;
    REQUIRE(qioAccessor.getHostRO()[0] == doctest::Approx(73.9465).epsilon(prec));
}
TEST_SUITE_END();

} // namespace Nextsim
