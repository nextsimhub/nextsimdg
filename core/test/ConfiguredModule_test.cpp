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

using HelpMap = int;
using Config = Nextsim::ConfiguredModule;

#include "include/Module.hpp"

#include <istream>
#include <memory>

namespace Test {
class ITest {};
class DefaultImpl : public ITest {};
class Impl1 : public ITest {};
class Impl2 : public ITest {};

class ITest2 {};
class AnotherImpl : public ITest2 {};
class DefaultImpl_2 : public ITest2 {};

}

// Module namespace to emulate how the names are defined with real Modules
namespace Module {
const std::string DEFAULTIMPL = "DefaultImpl";
const std::string IMPL1 = "Impl1";
const std::string IMPL2 = "Impl2";
const std::string ANOTHERIMPL = "AnotherImpl";

template <>
const Module<Test::ITest>::map& Module<Test::ITest>::functionMap()
{
    static const map theMap = {
        { DEFAULTIMPL, newImpl<Test::ITest, Test::DefaultImpl> },
        { IMPL1, newImpl<Test::ITest, Test::Impl1> },
        { IMPL2, newImpl<Test::ITest, Test::Impl2> },
    };
    return theMap;
}

template <>
Module<Test::ITest>::fn& Module<Test::ITest>::getGenerationFunction()
{
    static fn ptr = functionMap().at(DEFAULTIMPL);
    return ptr;
}

template <> std::string Module<Test::ITest>::moduleName() { return "ITest"; }

template <> HelpMap& Module<Test::ITest>::getHelpRecursive(HelpMap& map, bool getAll)
{
    return map;
}

template <>
const Module<Test::ITest2>::map& Module<Test::ITest2>::functionMap()
{
    static const map theMap = {
        { DEFAULTIMPL, newImpl<Test::ITest2, Test::DefaultImpl_2> },
        { ANOTHERIMPL, newImpl<Test::ITest2, Test::AnotherImpl> },
    };
    return theMap;
}

template <>
Module<Test::ITest2>::fn& Module<Test::ITest2>::getGenerationFunction()
{
    static fn ptr = functionMap().at(DEFAULTIMPL);
    return ptr;
}

template <> std::string Module<Test::ITest2>::moduleName() { return "ITest2"; }

template <> HelpMap& Module<Test::ITest2>::getHelpRecursive(HelpMap& map, bool getAll)
{
    return map;
}

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
    Module::Module<Test::ITest>::setImplementation(ConfiguredModule::getImpl("ITest"));
    REQUIRE(ConfiguredModule::getImpl("ITest2") == Module::ANOTHERIMPL);
    Module::Module<Test::ITest2>::setImplementation(ConfiguredModule::getImpl("ITest2"));

    // Tests Module::ImplementationNames, too
    ConfigMap cfgMap = ConfiguredModule::getAllModuleConfigurations();
    REQUIRE(cfgMap.size() == 2);
    REQUIRE(std::get<std::string>(cfgMap.at("ITest")) == Module::IMPL2);
    REQUIRE(std::get<std::string>(cfgMap.at("ITest2")) == Module::ANOTHERIMPL);
}

TEST_SUITE_END();

}
