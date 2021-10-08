/*!
 * @file CommandLineParser_test.cpp
 *
 * @date Oct 8, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "CommandLineParser.hpp"

#include <vector>
#include <string>

#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

namespace Nextsim {

// A class to work around the limitations on casting string literals to char*.
class ArgV {
public:
    ArgV(std::vector<std::string> vs)
    {
        nStrings = vs.size();
        ppc = new char*[nStrings + 1]; // Add 1 for the conventional final \0
        for (int i = 0; i < nStrings; ++i) {
            int len = vs[i].size();
            char* pc = new char[len + 1];
            for (int j = 0; j < len; ++j) {
                pc[j] = vs[i][j];
            }
            pc[len] = '\0';
            ppc[i] = pc;
        }
        char* pnull = new char[1];
        pnull[0] = '\0';
        ppc[nStrings] = pnull;
    }
    ~ArgV()
    {
        for (int i = 0; i < nStrings; ++i) {
            delete[] ppc[i];
        }
        delete[] ppc;
    }

    char** operator()()
    {
        return ppc;
    }

    int argc()
    {
        return nStrings + 1;
    }
private:
    int nStrings;
    char** ppc;
};

TEST_CASE("Parse config file names", "[CommandLineParser]")
{
    // Parse one file
    ArgV argv1({"nextsimdg", "--config-file", "config.cfg"});

    CommandLineParser clp1(argv1.argc(), argv1());
    std::vector<std::string> cfgs = clp1.getConfigFileNames();

    REQUIRE(cfgs.size() == 1);
    REQUIRE(cfgs[0] == std::string(argv1()[argv1.argc() - 2]));

}

} /* namespace Nextsim */
