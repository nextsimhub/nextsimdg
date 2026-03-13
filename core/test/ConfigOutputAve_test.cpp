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
    Nextsim::TwoDField cice(Nextsim::ModelArray::Type::TWOD);
    cice.resize();
    cice = 1.0;
    Nextsim::ModelState state;
    state.data["cice"] = cice;
    Nextsim::ConfigOutput co;
    co.outputState(state);
    auto oState = StructureFactory::getState();
    std::cout << oState.data.size() << std::endl;
}
TEST_SUITE_END();
