/*
 * @file ConfiguredModule_test.cpp
 *
 * @date Oct 21, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "ConfiguredModule.hpp"

#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include "ArgV.hpp"
#include "testClasses.hpp"
#include "Configurator.hpp"
#include "ModuleLoader.hpp"

namespace Nextsim {
TEST_CASE("Configure a module", "[Configurator, ModuleLoader]")
{
    // Create the fake command line, selecting the SNU2 albedo implementation
    ArgV argvee({"cmtest", "--ITest=Impl1"});

    Configurator::setCommandLine(argvee.argc(), argvee());

    Nextsim::ConfiguredModule::parseConfigurator();

    ITest& impler = ModuleLoader::getLoader().getImplementation<ITest>();

    REQUIRE(impler() == 1);
}
}
