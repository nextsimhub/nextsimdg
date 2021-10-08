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
        ppc = new char*[nStrings];
        for (int i = 0; i < nStrings; ++i) {
            int len = vs[i].size();
            char* pc = new char[len + 1];
            for (int j = 0; j < len; ++j) {
                pc[j] = vs[i][j];
            }
            pc[len] = '\0';
            ppc[i] = pc;
        }

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
private:
    int nStrings;
    char** ppc;
};

TEST_CASE("Parse config file names", "[CommandLineParser]")
{
    // Parse one file
    const int argc1 = 4;
    ArgV argv1({"nextsimdg", "--config-file", "config.cfg", ""});

    CommandLineParser clp1(argc1, argv1());
    std::vector<std::string> cfgs = clp1.getConfigFileNames();

    REQUIRE(cfgs.size() == 1);
    REQUIRE(cfgs[0] == std::string(argv1()[argc1 - 2]));
}

} /* namespace Nextsim */
