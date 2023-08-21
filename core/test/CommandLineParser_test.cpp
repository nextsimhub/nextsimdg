/*!
 * @file CommandLineParser_test.cpp
 *
 * @date Oct 8, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "ArgV.hpp"
#include "CommandLineParser.hpp"

#include <string>
#include <vector>

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

namespace Nextsim {

TEST_SUITE_BEGIN("CommandLineParser");
TEST_CASE("Parse config file names")
{
    // Parse one file
    ArgV argv1({ "nextsimdg", "--config-file", "config.cfg" });

    CommandLineParser clp1(argv1.argc(), argv1());
    std::vector<std::string> cfgs = clp1.getConfigFileNames();

    REQUIRE(cfgs.size() == 1);
    REQUIRE(cfgs[0] == std::string(argv1()[argv1.argc() - 1]));

    cfgs.clear();
    std::string finalFileName = "final.cfg";
    ArgV argv2({ "nextsimdg", "--config-file", "config.cfg", "--config-files", "test.cfg",
        "more.cfg", finalFileName });

    CommandLineParser clp2(argv2.argc(), argv2());
    cfgs = clp2.getConfigFileNames();

    REQUIRE(cfgs.size() == 4);
    REQUIRE(cfgs[cfgs.size() - 1] == finalFileName);
}
TEST_SUITE_END();

} /* namespace Nextsim */
