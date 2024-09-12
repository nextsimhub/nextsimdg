/*
 * @file ConfiguredModule_test.cpp
 *
 * @date Oct 21, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "ConfiguredModule.hpp"

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "ArgV.hpp"
#include "include/Configurator.hpp"
#include "include/Module.hpp"

#include <istream>
#include <memory>

// Module namespace to emulate how the names are defined with real Modules
namespace Module {
const std::string IMPL1 = "Impl1";
const std::string IMPL2 = "Impl2";
const std::string ANOTHERIMPL = "AnotherImpl";
}

namespace Nextsim {

TEST_SUITE_BEGIN("ConfiguredModule");
TEST_CASE("Configure a module")
{
    Configurator::clear();

    // Create the fake command line, selecting the Impl1 implementation
    ArgV argvee({"cmtest", "--Modules.ITest=Impl1"});

    Configurator::setCommandLine(argvee.argc(), argvee());

    REQUIRE(ConfiguredModule::getImpl("ITest") == Module::IMPL1);
}

TEST_CASE("Configure a module from a stream")
{
    Configurator::clear();
    std::stringstream config;
    config << "[Modules]" << std::endl
            << "ITest = Impl2" << std::endl
            << "ITest2 = AnotherImpl" << std::endl;
    std::unique_ptr<std::istream> pcstream(new std::stringstream(config.str()));
    Configurator::addStream(std::move(pcstream));

    REQUIRE(ConfiguredModule::getImpl("ITest") == Module::IMPL2);
    // Since no real Modules are instantiated, set the module names manually
    Module::ImplementationNames::set("ITest", Module::IMPL2);
    REQUIRE(ConfiguredModule::getImpl("ITest2") == Module::ANOTHERIMPL);
    Module::ImplementationNames::set("ITest2", Module::ANOTHERIMPL);

    ConfigMap cfgMap = ConfiguredModule::getAllModuleConfigurations();
    REQUIRE(cfgMap.size() == 2);
    REQUIRE(std::get<std::string>(cfgMap.at("ITest")) == Module::IMPL2);
    REQUIRE(std::get<std::string>(cfgMap.at("ITest2")) == Module::ANOTHERIMPL);
}



TEST_SUITE_END();

}
