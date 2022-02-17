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
#include "include/Configurator.hpp"
#include "include/Module.hpp"

#include <istream>
#include <memory>

// Module classes for the test classes
namespace Module {
class ITestModule : public Module<ITest>
{
};

const std::string IMPL1 = "Impl1";
const std::string IMPL2 = "Impl2";

template <>
Module<ITest>::map Module<ITest>::functionMap = {
        {IMPL1, newImpl<ITest, Impl1>},
        {IMPL2, newImpl<ITest, Impl2>},
};
template <> Module<ITest>::fn Module<ITest>::spf = functionMap.at(IMPL1);
template <> std::unique_ptr<ITest> Module<ITest>::staticInstance = std::move(Module<ITest>::spf());
template <> std::string Module<ITest>::moduleName() { return "ITest"; };
template <> ITest& getImplementation<ITest>() { return getImplTemplate<ITest, ITestModule>(); };
template <> void setImplementation<ITest>(const std::string& implName)
{
    setImplTemplate<ITestModule>(implName);
};
}

namespace Nextsim {
TEST_CASE("Configure a module", "[Configurator]")
{
    Configurator::clear();

    // Create the fake command line, selecting the Impl1 implementation
    ArgV argvee({"cmtest", "--Modules.ITest=Impl1"});

    Configurator::setCommandLine(argvee.argc(), argvee());

    ConfiguredModule::setConfiguredModules({
        {Module::Module<ITest>::moduleName(), Module::setImplementation<ITest>},
    });
    ConfiguredModule::parseConfigurator();

    ITest& impler = Module::getImplementation<ITest>();

    REQUIRE(impler() == 1);
}

TEST_CASE("Configure a module from a stream", "[Configurator]")
{
    Configurator::clear();
    std::stringstream config;
    config << "[Modules]" << std::endl
            << "ITest = Impl2" << std::endl;

    std::unique_ptr<std::istream> pcstream(new std::stringstream(config.str()));
    Configurator::addStream(std::move(pcstream));

    ConfiguredModule::setConfiguredModules({
        {Module::Module<ITest>::moduleName(), Module::setImplementation<ITest>},
    });
    ConfiguredModule::parseConfigurator();

    ITest& impler = Module::getImplementation<ITest>();
    REQUIRE(impler() == 2);
}

TEST_CASE("Don't configure a module from a stream", "[Configurator]")
{
    Configurator::clear();
    std::stringstream config;
    config << "[Modules]" << std::endl
            << "ITestNotReally = NotImpl2" << std::endl;

//    ModuleLoader::getLoader().setImplementation("ITest", "Impl2");
    Module::setImplementation<ITest>("Impl2");

    std::unique_ptr<std::istream> pcstream(new std::stringstream(config.str()));
    Configurator::addStream(std::move(pcstream));

    ConfiguredModule::setConfiguredModules({
        {Module::Module<ITest>::moduleName(), Module::setImplementation<ITest>},
    });
    // Parse the available modules. This should not change the implementation
    // to the default.
    ConfiguredModule::parseConfigurator();

    ITest& impler = Module::getImplementation<ITest>();
    REQUIRE(impler() == 2);
}

TEST_CASE("Configure a module with an incorrect name", "[Configurator]")
{
    Configurator::clear();
    std::stringstream config;
    config << "[Modules]" << std::endl
            << "ITest = Optometry" << std::endl;

    std::unique_ptr<std::istream> pcstream(new std::stringstream(config.str()));
    Configurator::addStream(std::move(pcstream));

    ConfiguredModule::setConfiguredModules({
        {Module::Module<ITest>::moduleName(), Module::setImplementation<ITest>},
    });
    // Should throw a domain_error as "Optometry" is not a valid implementation.
    REQUIRE_THROWS(ConfiguredModule::parseConfigurator());

}

}
