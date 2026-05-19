/*!
 * @author  Tim Spain <timothy.spain@nersc.no>
 */

#include <doctest/doctest.h>
#include <sstream>

#include "include/IFluxCalculation.hpp"
#include "include/ThermoIce0.hpp"

#include "include/Configurator.hpp"
#include "include/ConfiguredModule.hpp"
#include "include/IFreezingPoint.hpp"
#include "include/ModelArray.hpp"
#include "include/ModelComponent.hpp"
#include "include/NextsimModule.hpp"
#include "include/Time.hpp"

namespace Nextsim {

TEST_SUITE_BEGIN("ThermoIce0");
/*
 * Test that ice below the minimum ice threshold is eliminated.
 */
TEST_CASE("Threshold ice")
{
    ModelArray::setDimensions(ModelArray::Type::H, { 1, 1 });

    // Class derived from ModelComponent providing the physical data for the test
    class IceTemperatureData : public ModelComponent {
    public:
        IceTemperatureData()
            : hiceAccessor(getStore(), RW)
            , ciceAccessor(getStore(), RW)
            , hsnowAccessor(getStore(), RW)
            , sstAccessor(getStore(), RO)
            , sssAccessor(getStore(), RO)
            , tfAccessor(getStore(), RO)
            , snowAccessor(getStore(), RO)
            , mlbhcAccessor(getStore(), RO)
            , qioAccessor(getStore(), RW)
            , qowAccessor(getStore(), RW)
            , qiaAccessor(getStore(), RW)
            , dqia_dtAccessor(getStore(), RW)
            , sublAccessor(getStore(), RW)
            , penSWAccessor(getStore(), RW)
        {
        }
        std::string getName() const override { return "IceTemperatureData"; }

        void setData(const ModelState::DataMap&) override
        {
            ciceAccessor.getHostRW()[0] = 0.99;
            HField& hice = hiceAccessor.getHostRW();
            hice[0] = 0.001;
            hice0[0] = hice[0];
            hsnowAccessor.getHostRW()[0] = 0.;
            HField& sss = sssAccessor.getHostRW();
            sss[0] = 32.;
            sstAccessor.getHostRW()[0] = Module::getImplementation<IFreezingPoint>()(sss[0]);
            snowAccessor.getHostRW()[0] = 0.;
            tfAccessor.getHostRW()[0] = Module::getImplementation<IFreezingPoint>()(sss[0]);
            mlbhcAccessor.getHostRW()[0] = 4.29151e7;
            qioAccessor.getHostRW()[0] = 0.;
            qowAccessor.getHostRW()[0] = 0;
            qiaAccessor.getHostRW()[0] = 0;
            dqia_dtAccessor.getHostRW()[0] = 0;
            sublAccessor.getHostRW()[0] = 0;
            penSWAccessor.getHostRW()[0] = 0;
        }

        ModelArrayAccessor<Shared::H_ICE_DG, RW> hiceAccessor;
        HField hice0;
        ModelArrayAccessor<Shared::C_ICE_DG, RW> ciceAccessor;
        ModelArrayAccessor<Shared::H_SNOW_DG, RW> hsnowAccessor;
        ModelArrayAccessor<Protected::SST, RW> sstAccessor;
        ModelArrayAccessor<Protected::SSS, RW> sssAccessor;
        ModelArrayAccessor<Protected::TF, RW> tfAccessor;
        ModelArrayAccessor<Protected::SNOW, RW> snowAccessor;
        ModelArrayAccessor<Protected::ML_BULK_CP, RW>
            mlbhcAccessor; // Mixed layer bulk heat capacity
        ModelArrayAccessor<Shared::Q_IO, RW> qioAccessor;
        ModelArrayAccessor<Shared::Q_OW, RW> qowAccessor;
        ModelArrayAccessor<Shared::Q_IA, RW> qiaAccessor;
        ModelArrayAccessor<Shared::DQIA_DT, RW> dqia_dtAccessor;
        ModelArrayAccessor<Shared::SUBLIM, RW> sublAccessor;
        ModelArrayAccessor<Shared::Q_PEN_SW, RW> penSWAccessor;
    } atmoState;
    atmoState.setData(ModelState::DataMap());

    // Supply the atmospheric boundary arrays, without an entire
    // IAtmosphereBoundary implementation
    ModelArrayAccessor<Shared::Q_IC, RW> qicAccessor(ModelComponent::getStore(), RW);
    HField& qic = qicAccessor.getHostRW();
    qic.reinitialize();

    ModelArrayAccessor<Shared::Q_SW_BASE, RW> qswbaseAccessor(ModelComponent::getStore(), RW);
    HField& qswbase = qswbaseAccessor.getHostRW();
    qswbase.reinitialize();

    // An implementation of IFluxCalculation that returns zero fluxes
    class FluxData : public IFluxCalculation {
    public:
        FluxData()
            : IFluxCalculation()
        {
        }
        std::string getName() const override { return "FluxData"; }

        void setData(const ModelState::DataMap&) override
        {
            qowAccessor.getHostRW()[0] = 0;
            qiaAccessor.getHostRW()[0] = 0;
            dqia_dtAccessor.getHostRW()[0] = 0;
            sublAccessor.getHostRW()[0] = 0;
            penSWAccessor.getHostRW()[0] = 0;
        }

        void update(const TimestepTime&) override { }
    } fluxData;

    fluxData.setData(ModelState::DataMap());

    TimestepTime tst = { TimePoint("2000-01-01T00:00:00"), Duration(600) };
    ThermoIce0 ti0t;
    ti0t.configure();

    HField tSurf;
    tSurf = -9.;
    ti0t.setData({ { tsurfName, tSurf } });
    ti0t.update(tst);

    ModelArrayAccessor<Shared::H_ICE_DG> hiceAccessor(ModelComponent::getStore());
    const HField& hice = hiceAccessor.getHostRO();

    // So little ice should be reduced to zero
    REQUIRE(hice[0] == 0.);
    ModelArrayAccessor<Shared::C_ICE_DG> ciceAccessor(ModelComponent::getStore());
    const HField& cice = ciceAccessor.getHostRO();

    REQUIRE(cice[0] == 0.);
}
TEST_SUITE_END();
}
