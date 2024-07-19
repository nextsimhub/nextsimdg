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
#include "include/IStructure.hpp"
#include "include/ModelComponent.hpp"
#include "include/Module.hpp"
#include "include/StructureFactory.hpp"
#include "include/UnescoFreezing.hpp"
#include "include/constants.hpp"

#include <filesystem>
#include <sstream>

extern template class Module::Module<Nextsim::IOceanBoundary>;

namespace Nextsim {

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
#define isDG 1 // remove me
#if isDG
TEST_CASE("PrognosticData write test, including DG components")
{
    std::string filename = "dg.nc";

    // Set the structure type
    Module::setImplementation<IStructure>("Nextsim::ParametricGrid");
    // Set the dimensions using the newer per-dimension interface.
    ModelArray::setDimension(ModelArray::Dimension::X, 2);
    ModelArray::setDimension(ModelArray::Dimension::Y, 2);
    ModelArray::setDimension(ModelArray::Dimension::Z, 1);
    ModelArray::setDimension(ModelArray::Dimension::XVERTEX, ModelArray::size(ModelArray::Dimension::X) + 1);
    ModelArray::setDimension(ModelArray::Dimension::YVERTEX, ModelArray::size(ModelArray::Dimension::Y) + 1);

    REQUIRE(ModelArray::nComponents(ModelArray::Type::DG) > 1);

    // Set the data, fill out the arrays
    ModelArray baseData(ModelArray::Type::H);
    ModelArray zData(ModelArray::Type::Z);
    ModelArray dgData(ModelArray::Type::DG);
    ModelArray latitude(ModelArray::Type::H);
    ModelArray longitude(ModelArray::Type::H);
    ModelArray coords(ModelArray::Type::VERTEX);
    baseData.resize();
    zData.resize();
    dgData.resize();
    coords.resize();
    size_t xMul = 1000;
    size_t yMul = 100;
    size_t resn = 2;
    for (size_t j = 0; j < ModelArray::size(ModelArray::Dimension::Y); ++j) {
        for (size_t i = 0; i < ModelArray::size(ModelArray::Dimension::X); ++i) {
            size_t c = i * xMul + j * yMul;
            size_t idx = baseData.indexFromLocation({i, j});
            baseData[idx] = c;
            zData.zIndexAndLayer(idx, 0) = c;
            latitude(idx) = resn * j + 0.5 * resn;
            longitude(idx) = resn * i + 0.5 * resn;
            switch(ModelArray::nComponents(ModelArray::Type::DG)) {
            case(3):
                    dgData.components(idx) << c, c + 1, c + 2;
            break;
            case(6):
                    dgData.components(idx) << c, c + 1, c + 2, c + 3, c + 4, c + 5;
            break;
            default:
                REQUIRE_MESSAGE(false, "Unimplemented number of DG components");
            }
        }
    }

    // Loop for coords
    for (size_t j = 0; j < ModelArray::size(ModelArray::Dimension::YVERTEX); ++j) {
        for (size_t i = 0; i < ModelArray::size(ModelArray::Dimension::XVERTEX); ++i) {
            size_t c = i * xMul + j * yMul;
            size_t idx = coords.indexFromLocation({i, j});
            coords.components(idx) << resn * i, resn * j;
        }
    }

    // Create the data map
    std::map<std::string, size_t> offsets = { {ciceName, 10}, {hiceName, 20}, {hsnowName, 30}, {sssName, 40}, {sstName, 50}, {ticeName, 60}, {uName, 70}, {vName, 80}};
    ModelState::DataMap inData = {
            { ciceName, dgData + offsets[ciceName] },
            { hiceName, dgData + offsets[hiceName] },
            { hsnowName, baseData + offsets[hsnowName] },
            { sssName, baseData + offsets[sssName] },
            { sstName, baseData + offsets[sstName] },
            { ticeName, zData + offsets[ticeName] },
            { uName, baseData + offsets[uName] },
            { vName, baseData + offsets[vName] },
            // grid coordinates and mask
            { maskName, baseData * 0 + 1},
            { gridAzimuthName, baseData * 0 },
            { latitudeName, latitude },
            { longitudeName, longitude },
            { coordsName, coords },
    };

    // Create the writing classes, data and metadata
    ModelState state = {inData, {} };
    PrognosticData pData;
    pData.configure();
    pData.setData(inData);
    ModelMetadata meta;
    meta.setTime(TimePoint("2010-01-01T00:00:00Z"));
    meta.extractCoordinates(state);
    // Do the write
    pData.writeRestartFile(filename, meta);

    // Read the data back
    ModelState readState = StructureFactory::stateFromFile(filename);
    // Test the contents of the data
    ModelState::DataMap& readData = readState.data;
    // Check for the coordinates
    REQUIRE(readData.count(coordsName) > 0);
    REQUIRE(readData.count(latitudeName) > 0);
    // Check the type and contents of cice
    REQUIRE(readData.count(ciceName) > 0);
    ModelArray& cice = readData.at(ciceName);
    REQUIRE(cice.getType() == ModelArray::Type::DG);
    size_t testComponent = 4;
    REQUIRE(cice.components(cice.indexFromLocation({1, 1}))[testComponent] == 1 * xMul + 1 * yMul + testComponent + offsets[ciceName]);

    std::filesystem::remove(filename);


}

#endif // isDG
TEST_SUITE_END();

} /* namespace Nextsim */
