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

    class NewIceData : public ModelComponent {
    public:
        NewIceData()
        {
            registerProtectedArray(ProtectedArray::H_ICE, &hice_old);
            registerProtectedArray(ProtectedArray::SST, &sst);
            registerProtectedArray(ProtectedArray::SSS, &sss);
            registerProtectedArray(ProtectedArray::TF, &tf);
            registerProtectedArray(ProtectedArray::SNOW, &snow);
            registerProtectedArray(ProtectedArray::ML_BULK_CP, &mlbhc);

            registerSharedArray(SharedArray::H_ICE, &hice);
            registerSharedArray(SharedArray::C_ICE, &cice);
            registerSharedArray(SharedArray::H_SNOW, &hsnow);
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
            hice_old[0] = 0.2;
            sst[0] = -1.2;
            sss[0] = 32.;
            snow[0] = 0.;
            tf[0] = -1.8;

            hice[0] = hice_old[0];
            cice[0] = 0.5;
            hsnow[0] = 0;
            qow[0] = 124.689;
            qio[0] = 124.689;
            qic[0] = 0.; // Get the data value
            qia[0] = 305.288;
            dqia_dt[0] = 4.5036;
            subl[0] = 0;

            tice[0] = -2;
        }

        HField hice_old;
        HField sst;
        HField sss;
        HField tf;
        HField snow;
        HField mlbhc; // Mixed layer bulk heat capacity

        HField hice;
        HField cice;
        HField hsnow;
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

    REQUIRE(newice[0] == 0.0258264);
}

}
