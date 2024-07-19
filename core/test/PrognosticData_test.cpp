/*!
 * @file PrognosticData_test.cpp
 *
 * @date 7 Sep 2023
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/PrognosticData.hpp"

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "include/PrognosticData.hpp"

#include "include/ConfiguredModule.hpp"
#include "include/ModelComponent.hpp"
#include "include/Module.hpp"
#include "include/UnescoFreezing.hpp"
#include "include/constants.hpp"

#include <sstream>

extern template class Module::Module<Nextsim::IOceanBoundary>;

namespace Nextsim {

void PrognosticData::writeRestartFile(const std::string& filePath, const ModelMetadata&) const
{
}

TEST_SUITE_BEGIN("PrognosticData");
TEST_CASE("PrognosticData call order test")
{
    ModelArray::setDimensions(ModelArray::Type::H, { 1, 1 });
    ModelArray::setDimensions(ModelArray::Type::Z, { 1, 1, 1 });

    std::stringstream config;
    config << "[Modules]" << std::endl;
    config << "AtmosphereBoundaryModule = Nextsim::ConfiguredAtmosphere" << std::endl;
    config << std::endl;
    config << "[ConfiguredAtmosphere]" << std::endl;
    config << "t_air = 3" << std::endl;
    config << "t_dew = 2" << std::endl;
    config << "pmsl = 100000" << std::endl;
    config << "sw_in = 50" << std::endl;
    config << "lw_in = 330" << std::endl;
    config << "snow = 0" << std::endl;
    config << "rainfall = 0" << std::endl;
    config << "wind_speed = 5" << std::endl;

    std::unique_ptr<std::istream> pcstream(new std::stringstream(config.str()));
    Configurator::addStream(std::move(pcstream));

    ConfiguredModule::parseConfigurator();

    class OceanData : public IOceanBoundary {
    public:
        OceanData()
            : IOceanBoundary()
        {
        }
        void setData(const ModelState::DataMap& state) override { IOceanBoundary::setData(state); }
        void updateBefore(const TimestepTime& tst) override
        {
            UnescoFreezing uf;
            sst = -1.;
            sss = 32.;
            mld = 10.25;
            tf = uf(sss[0]);
            cpml = Water::cp * Water::rho * mld[0];
            u = 0;
            v = 0;
        }
        void updateAfter(const TimestepTime& tst) override { }
    } ocnBdy;
    ocnBdy.setData(ModelState().data);

    Module::Module<IOceanBoundary>::setExternalImplementation(
        Module::newImpl<IOceanBoundary, OceanData>);

    HField zeroData;
    zeroData.resize();
    zeroData[0] = 0.;
    ZField zeroDataZ;
    zeroDataZ.resize();
    zeroDataZ[0] = 0.;

    ModelState::DataMap initialData = {
        { "cice", zeroData },
        { "hice", zeroData },
        { "hsnow", zeroData },
        { "tice", zeroDataZ },
    };

    PrognosticData pData;
    pData.configure();
    pData.setData(initialData);
    TimestepTime tst = { TimePoint("2000-01-01T00:00:00Z"), Duration("P0-0T0:10:0") };
    pData.update(tst);

    ModelArrayRef<Shared::Q_OW> qow(ModelComponent::getStore());

    double prec = 1e-5;
    // Correct value
    REQUIRE(qow[0] == doctest::Approx(-109.923).epsilon(prec));
    // Value if pAtmBdy->update and pOcnBdy->updateBefore are switched in PrognosticData::update
    REQUIRE(qow[0] != doctest::Approx(-92.1569).epsilon(prec));
}
TEST_SUITE_END();

} /* namespace Nextsim */
