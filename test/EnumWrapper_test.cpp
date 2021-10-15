/*!
 * @file EnumWrapper_test.cpp
 *
 * @date Oct 15, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "EnumWrapper.hpp"
#include "ArgV.hpp"

#include <string>
#include <sstream>
#include <boost/program_options.hpp>

#define CATCH_CONFIG_MAIN
//#include <catch2/catch.hpp>
#include "/opt/home/include/catch2/catch.hpp"

namespace Nextsim {

enum class Colour {
    red,
    green,
    blue,
};

TEST_CASE("Configure a wrapped enum", "[EnumWrapper]")
{
    std::string targetColourString = "red";

    EnumWrap::EnumWrapper<Colour>::setMap({
        {targetColourString, Colour::red},
        {"green", Colour::green},
        {"blue", Colour::blue},
    });

    std::string colourName = "option.colour";
    Colour targetColour = Colour::red;

    std::stringstream cfgText;
    cfgText << "[option]" << std::endl
            << "colour = " << targetColourString << std::endl;

    boost::program_options::options_description opt("EnumWrapper test options");
    opt.add_options()
            (colourName.c_str(), boost::program_options::value<EnumWrap::EnumWrapper<Colour>>(), "What is you favourite colour?")
            ;
    boost::program_options::variables_map vm;
    boost::program_options::store(boost::program_options::parse_config_file(cfgText, opt, true), vm);

    Colour col = vm[colourName].as<EnumWrap::EnumWrapper<Colour>>();
    REQUIRE(col == targetColour);
}
}
