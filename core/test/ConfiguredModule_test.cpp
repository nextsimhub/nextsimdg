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
#include "moduleTestClasses.hpp"
#include "Configurator.hpp"
#include "ModuleLoader.hpp"

#include <istream>
#include <memory>

namespace Nextsim {
TEST_CASE("Configure a module", "[Configurator, ModuleLoader]")
{
    Configurator::clear();

    // Create the fake command line, selecting the Impl1 implementation
    ArgV argvee({"cmtest", "--Modules.ITest=Impl1"});

    Configurator::setCommandLine(argvee.argc(), argvee());

    ConfiguredModule::parseConfigurator();

    ITest& impler = ModuleLoader::getLoader().getImplementation<ITest>();

    REQUIRE(impler() == 1);
}

TEST_CASE("Configure a module from a stream", "[Configurator, ModuleLoader]")
{
    Configurator::clear();
    std::stringstream config;
    config << "[Modules]" << std::endl
            << "ITest = Impl2" << std::endl;

    std::unique_ptr<std::istream> pcstream(new std::stringstream(config.str()));
    Configurator::addStream(std::move(pcstream));

    ConfiguredModule::parseConfigurator();

    ITest& impler = ModuleLoader::getLoader().getImplementation<ITest>();
    REQUIRE(impler() == 2);
}


}
