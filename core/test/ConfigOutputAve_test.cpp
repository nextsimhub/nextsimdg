/*!
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include <string>
#include <iostream>

/**********************************************************************
 *
 * Mock classes
 *
 *********************************************************************/
#include "include/StructureFactory.hpp"

#include "modules/DiagnosticOutputModule/include/ConfigOutput.hpp"

const std::map<std::string, std::string> IDiagnosticOutput::externalNames = {
};

/**********************************************************************
 *
 * Test code
 *
 *********************************************************************/
const auto nx = 1;
const auto ny = 1;

TEST_SUITE_BEGIN("ConfigOutput");
TEST_CASE("Test single output") {
    Nextsim::ModelArray::setDimensions(Nextsim::ModelArray::Type::TWOD, {nx, ny});
    Nextsim::TwoDField cice(Nextsim::ModelArray::Type::TWOD);
    // Clear the mock StructureFactory class
    StructureFactory::getState().data.clear();
    cice.resize();
    Nextsim::ModelState state;
    Nextsim::ConfigOutput co;
    co.configure();
    // Set the time in the mock ModelMetadata class
    ModelMetadata::getInstance().startTime() = Nextsim::TimePoint("2000-01-01T00:00:00Z");
    ModelMetadata::getInstance().time() = Nextsim::TimePoint("2000-01-01T01:00:00Z");
    const std::string cice_name = "cice";
    double cice_val = 1.;
    cice = cice_val;
    state.data[cice_name] = cice;
    co.outputState(state);
    auto oState = StructureFactory::getState();
    REQUIRE(oState.data.size() == 1);
    REQUIRE(oState.data[cice_name][0] == cice_val);
}

TEST_CASE("Averaging output") {
    Nextsim::ModelArray::setDimensions(Nextsim::ModelArray::Type::TWOD, {nx, ny});
    Nextsim::TwoDField cice(Nextsim::ModelArray::Type::TWOD);
    StructureFactory::getState().data.clear();
    cice.resize();
    Nextsim::ModelState state;
    Nextsim::ConfigOutput co;
    co.configure();
    ModelMetadata::getInstance().startTime() = Nextsim::TimePoint("2000-01-01T00:00:00Z");
    const std::string cice_name = "cice";
    double accum = 0.;
    const size_t nt = 4;
    for (auto i = 1; i <= nt; ++i) {
        auto t = Nextsim::TimePoint(std::string("2000-01-01T0")+std::to_string((i*15)/60)+":"+std::to_string((i*15) % 60)+":00Z");
        // Set the step time in the mock ModelMetadata class
        ModelMetadata::getInstance().time() = t;
        double cice_val = i;
        cice = cice_val;
        accum += cice_val;
        state.data[cice_name] = cice;
        co.outputState(state);

    }
    accum /= nt;
    auto oState = StructureFactory::getState();
    REQUIRE(oState.data.size() == 1);
    REQUIRE(oState.data[cice_name][0] == accum);
}

TEST_CASE("Model start before output start") {
    Nextsim::ModelArray::setDimensions(Nextsim::ModelArray::Type::TWOD, {nx, ny});
    Nextsim::TwoDField cice(Nextsim::ModelArray::Type::TWOD);
    StructureFactory::getState().data.clear();
    cice.resize();
    Nextsim::ModelState state;
    Nextsim::ConfigOutput co;
    co.configure();
    ModelMetadata::getInstance().startTime() = Nextsim::TimePoint("1999-12-31T00:00:00Z");
    const std::string cice_name = "cice";
    auto oState = StructureFactory::getState();
    REQUIRE(oState.data.size() == 0);
    // Initial steps
    size_t nt = 96;
    for (auto i = 1; i <= nt; ++i) {
        std::string day_str = (i < 96) ? "1999-12-31" : "2000-01-01";
        auto t = Nextsim::TimePoint(std::string(day_str+"T0")+std::to_string((i*15)/60)+":"+std::to_string((i*15) % 60)+":00Z");
        // Set the step time in the mock ModelMetadata class
        ModelMetadata::getInstance().time() = t;
        double cice_val = i;
        cice = cice_val;
        state.data[cice_name] = cice;
        co.outputState(state);
    }
    oState = StructureFactory::getState();
    REQUIRE(oState.data.size() == 0);

    double accum = 0.;
    size_t nt_off = 96; // Number of steps in the previous day
    nt = 4;
    for (auto i = 1; i <= nt; ++i) {
        auto t = Nextsim::TimePoint(std::string("2000-01-01T0")+std::to_string((i*15)/60)+":"+std::to_string((i*15) % 60)+":00Z");
        // Set the step time in the mock ModelMetadata class
        ModelMetadata::getInstance().time() = t;
        double cice_val = i + nt_off;
        cice = cice_val;
        accum += cice_val;
        state.data[cice_name] = cice;
        co.outputState(state);

    }
    accum /= nt;
    oState = StructureFactory::getState();
    REQUIRE(oState.data.size() == 1);
    REQUIRE(oState.data[cice_name][0] == accum);

}
TEST_SUITE_END();
