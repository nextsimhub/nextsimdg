/*!
 * @file ConfigOutputAve_test.cpp
 *
 * @date Mar 12, 2026
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
    cice.resize();
    Nextsim::ModelState state;
    Nextsim::ConfigOutput co;
    co.configure();
    auto t = Nextsim::TimePoint("2000-01-01T01:00:00Z");
    ModelMetadata::getInstance().time() = t;
    std::string cice_name = "cice";
    double cice_val = 1.;
    cice = cice_val;
    state.data[cice_name] = cice;
    co.outputState(state);
    auto oState = StructureFactory::getState();
    REQUIRE(oState.data.size() == 1);
    REQUIRE(oState.data[cice_name][0] == cice_val);
}
TEST_SUITE_END();
