/*
 * @file ModuleLoader_test.cpp
 *
 * @date Sep 23, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "ModuleLoader.hpp"

#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

# include "moduleTestClasses.hpp"
TEST_CASE("Basic module loading test", "[ModuleLoader]")
{
    ModuleLoader& ldr = ModuleLoader::getLoader();

    REQUIRE(ldr.listModules().size() == 1);
    REQUIRE(ldr.listImplementations("ITest").size() == 2);

    ldr.setImplementation("ITest", "Impl1");

    Impl1 i1;

    REQUIRE(typeid(i1) == typeid(*(ldr.getInstance<ITest>())));
}
