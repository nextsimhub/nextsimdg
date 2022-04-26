/*!
 * @file IceGrowth_test.cpp
 *
 * @date Apr 8, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>
#include <sstream>

#include "include/IceGrowth.hpp"

#include "include/Configurator.hpp"
#include "include/ConfiguredModule.hpp"
#include "include/IFreezingPoint.hpp"
#include "include/IFreezingPointModule.hpp"
#include "include/ModelArray.hpp"
#include "include/ModelArrayRef.hpp"
#include "include/ModelComponent.hpp"
#include "include/Time.hpp"

namespace Nextsim {

TEST_CASE("New ice formation", "[IceGrowth]")
{
    ModelArray::setDimensions(ModelArray::Type::H, { 1 });
    ModelArray::setDimensions(ModelArray::Type::Z, { 1, 1 });

    std::stringstream config;
    config << "[Modules]" << std::endl;
    config << "Nextsim::IFreezingPoint = Nextsim::UnescoFreezing" << std::endl;

    std::unique_ptr<std::istream> pcstream(new std::stringstream(config.str()));
    Configurator::addStream(std::move(pcstream));

    ConfiguredModule::parseConfigurator();

    class NewIceData : public ModelComponent {
    public:
        NewIceData()
        {
            registerProtectedArray(ProtectedArray::H_ICE, &hice);
            registerProtectedArray(ProtectedArray::C_ICE, &cice);
            registerProtectedArray(ProtectedArray::H_SNOW, &hsnow);
            registerProtectedArray(ProtectedArray::SST, &sst);
            registerProtectedArray(ProtectedArray::SSS, &sss);
            registerProtectedArray(ProtectedArray::TF, &tf);
            registerProtectedArray(ProtectedArray::SNOW, &snow);
            registerProtectedArray(ProtectedArray::ML_BULK_CP, &mlbhc);

            registerSharedArray(SharedArray::Q_OW, &qow);
            registerSharedArray(SharedArray::Q_IC, &qic);
            registerSharedArray(SharedArray::Q_IO, &qio);
            registerSharedArray(SharedArray::Q_IA, &qia);
            registerSharedArray(SharedArray::DQIA_DT, &dqia_dt);
            registerSharedArray(SharedArray::SUBLIM, &subl);

            registerSharedArray(SharedArray::T_ICE, &tice);
        }
        std::string getName() const override { return "NewIceData"; }

        void setData(const ModelState&) override
        {
            hice[0] = 0.2;
            cice[0] = 0.5;
            hsnow[0] = 0;
            sst[0] = -1.5;
            sss[0] = 32.;
            snow[0] = 0.;
            tf[0] = Module::getImplementation<IFreezingPoint>()(sss[0]);
            mlbhc[0] = 4.29151e7;

            qow[0] = 307.546;
            qio[0] = 124.689;
            qic[0] = 2.53124;
            qia[0] = 305.288;
            dqia_dt[0] = 4.5036;
            subl[0] = 0;

            tice[0] = -2;
        }

        HField hice;
        HField cice;
        HField hsnow;
        HField sst;
        HField sss;
        HField tf;
        HField snow;
        HField mlbhc; // Mixed layer bulk heat capacity

        HField qow;
        HField qic;
        HField qio;
        HField qia;
        HField dqia_dt;
        HField subl;

        ZField tice;
        ModelState getState() const override { return ModelState(); }
        ModelState getState(const OutputLevel&) const override { return getState(); }

    } iceData;

    iceData.setData(ModelState());

    TimestepTime tst = {0, 86400};
    IceGrowth ig;
    ig.configure();
    ig.update(tst);

    ModelArrayRef<ModelComponent::SharedArray::NEW_ICE, RO> newice;

    REQUIRE(newice[0] == Approx(0.0258264).epsilon(1e-5));
}

}
